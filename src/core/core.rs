use crate::core::helper::{
    average_vec_u16, is_all_ones, is_all_zeros, median, merge_u32_to_u64, percentile, to_string1,
    wrapper_stats, Feature,
};

use crate::alignment::pack::matrix_pack_wrapper;
use crate::r#mod::input_data::FileData;
use bitvec::prelude::*;
use gfa_reader::NCGfa;
use hashbrown::HashSet;
use packing_lib::convert::convert_helper::Method;
use packing_lib::convert::helper::median_vec_u16_16;
use std::fmt::Debug;
use std::fs::File;
use std::io::{BufWriter, Write};

#[derive(Debug, Clone, Eq, PartialEq)]
/// Core data structure
///
/// SNP-major representation of the genotypes (also in memory) - Reading data takes longer here
/// Can represent two haplotypes (diploid) or one haplotype (haploid)
pub struct MatrixWrapper {
    pub shape: (usize, usize),
    pub matrix_u16: Vec<Vec<u16>>,
    pub matrix_bit: Vec<BitVec<u8, Lsb0>>, //Vec<Vec<bool>>

    // Check if node, edges, dirnode, or alignment
    pub feature: Feature,     // Feature
    pub geno_names: Vec<u64>, // Name of all
    pub window_number: Vec<u32>,
    pub window_size: usize,                // Window number
    pub sample_names: Vec<String>,         // Sample names (same order as in the matrix)
    pub sample_index_u16: Vec<[usize; 2]>, // Sample index

    // Plink specific stuff
    pub fam_entries: Vec<String>, // Fam entries
    pub bim_entries: Vec<String>, // Bim entries
}

impl MatrixWrapper {
    /// Dummy initialization
    pub fn new() -> Self {
        Self {
            shape: (0, 0),
            matrix_u16: Vec::new(),
            matrix_bit: Vec::new(),
            feature: Feature::Node,

            // SNP
            geno_names: Vec::new(),
            window_number: Vec::new(),
            window_size: 0,
            bim_entries: Vec::new(),

            // Fam
            sample_names: Vec::new(),
            fam_entries: Vec::new(),
            sample_index_u16: Vec::new(),
        }
    }

    /// Initialize the SNP index
    pub fn make_index(&mut self, data: &NCGfa<()>, t: Feature) {
        //let mut bb = HashMap::with_hasher(BuildHasherDefault::<NoHashHasher<u32>>::default());
        let mut geno_names = Vec::new();
        match t {
            Feature::Node => {
                for (_i, x) in data.nodes.iter().enumerate() {
                    geno_names.push(x.id as u64);
                    //bb.insert(x.id as u64, i);
                }
            }
            Feature::DirNode => {
                if data.edges.is_some() {
                    let value = data.edges.as_ref().unwrap();
                    let mut edd = HashSet::new();
                    for x in value.iter() {
                        edd.insert(x.from as u64 * 2 + x.from_dir as u64);
                        edd.insert(x.to as u64 * 2 + x.to_dir as u64);
                    }
                    let mut edd2 = edd.into_iter().collect::<Vec<u64>>();
                    edd2.sort();

                    for (_i, x) in edd2.iter().enumerate() {
                        geno_names.push(*x);
                    }
                }
            }
            Feature::Edge => {
                if data.edges.is_some() {
                    let value = data.edges.as_ref().unwrap();
                    for (_i, x) in value.iter().enumerate() {
                        let (v1, v2, v3, v4) = (x.from, x.from_dir, x.to, x.to_dir);
                        let u1 = v1 * 2 + v2 as u32;
                        let u2 = v3 * 2 + v4 as u32;
                        let uu = merge_u32_to_u64(u1, u2);
                        geno_names.push(uu);
                    }
                }
            }
            _ => {}
        }
        // Sort it, otherwise does not work
        geno_names.sort();
        self.geno_names = geno_names;
    }

    pub fn make_thresh(&self, absolute: u32, relative: u16, method: Method) -> Vec<f32> {
        if absolute != 0 {
            return vec![absolute as f32; self.matrix_u16.len()];
        }
        match method {
            Method::Mean => {
                let mut aa = Vec::new();
                for x in self.matrix_u16.iter() {
                    aa.push(average_vec_u16(x, relative) as f32);
                }
                aa
            }
            Method::Median => {
                let mut aa = Vec::new();
                for x in self.matrix_u16.iter() {
                    aa.push(median(x, relative) as f32);
                }
                aa
            }
            Method::Percentile => {
                let mut aa = Vec::new();
                for x in self.matrix_u16.iter() {
                    aa.push(percentile(x, relative as f64) as f32);
                }
                aa
            }

            _ => {
                let mut aa = Vec::new();
                for x in self.matrix_u16.iter() {
                    aa.push(*x.iter().max().unwrap() as f32);
                }
                aa
            }
        }
    }

    pub fn make_bin_row(&mut self, relative: &Vec<f32>) {
        let mut matrix_bin = Vec::new();
        for (val, re) in self.matrix_u16.iter().zip(relative.iter()) {
            let mut biit = BitVec::with_capacity(val.len());
            for y in val.iter() {
                if *y as f64 >= *re as f64 {
                    biit.push(true);
                } else {
                    biit.push(false);
                }
            }
            matrix_bin.push(biit);
        }

        self.matrix_bit = matrix_bin;
    }

    /// Remove those entries which hold no information (all zero or all one)
    pub fn remove_non_info(&mut self) {
        let mut write_index = 0;
        for read_index in 0..self.geno_names.len() {
            if !is_all_ones(&self.matrix_bit[read_index])
                == !is_all_zeros(&self.matrix_bit[read_index])
            {
                // Retain elements that satisfy the condition
                if write_index != read_index {
                    // Move the elements to their new positions
                    self.matrix_bit.swap(write_index, read_index);
                    self.geno_names.swap(write_index, read_index);
                }

                write_index += 1;
            }
        }
        self.matrix_bit.truncate(write_index);
        self.geno_names.truncate(write_index);
    }

    /// Remove samples from the dataset
    pub fn remove_samples(&mut self, to_be_removed: &Vec<String>) {
        let mut remove_index: Vec<usize> = Vec::new();
        let samples_to_be_removed: HashSet<String> = to_be_removed.iter().cloned().collect();
        for (i, x) in self.sample_names.iter().enumerate() {
            if samples_to_be_removed.contains(x) {
                remove_index.push(i - remove_index.len());
            }
        }

        self.sample_names
            .retain(|x| !samples_to_be_removed.contains(x));
        if !self.fam_entries.is_empty() {
            self.fam_entries.retain(|x| {
                !samples_to_be_removed
                    .contains(&x.split_whitespace().collect::<Vec<&str>>()[0].to_string())
            });
        }
        assert_eq!(self.sample_names.len(), self.fam_entries.len());
        let mut i = 0;
        while i < self.matrix_bit.len() {
            let a = &mut self.matrix_bit[i];
            let mut f = 0;
            for x in remove_index.iter() {
                a.remove(*x * 2 - f * 2);
                a.remove(*x * 2 - f * 2 + 1);
                f += 1
            }
            i += 1;
        }
    }

    /// Remove entries from a file
    ///
    /// Remove:
    ///    - Entries
    pub fn remove_feature(&mut self, d: &FileData) {
        let mut write_index = 0;
        let mut i = 0;
        let mut j = 0;
        while i < self.geno_names.len() && j < d.data.len() {
            if self.geno_names[i] != d.data[j] {
                if write_index != i {
                    // Move the elements to their new positions
                    self.matrix_bit.swap(write_index, i);
                    self.geno_names.swap(write_index, i);
                    self.bim_entries.swap(write_index, i);
                    self.window_number.swap(write_index, i);
                }
                write_index += 1;
                if self.geno_names[i] < d.data[j] {
                    i += 1;
                } else {
                    j += 1;
                }
            } else {
                i += 1;
                j += 1;
            }
        }
        if i != self.geno_names.len() {
            for x in i..self.geno_names.len() {
                if write_index != x {
                    // Move the elements to their new positions
                    self.matrix_bit.swap(write_index, x);
                    self.geno_names.swap(write_index, x);
                    self.bim_entries.swap(write_index, x);
                    self.window_number.swap(write_index, x);
                }
                write_index += 1;
            }
        }

        println!("Write index: {:?}", write_index);
        self.matrix_bit.truncate(write_index);
        self.geno_names.truncate(write_index);
        self.bim_entries.truncate(write_index);
        self.window_number.truncate(write_index);
    }



    pub fn remove_feature_from_index(&mut self, d: &FileData) {
        let mut write_index = 0;
        let mut i = 0;
        let mut j = 0;
        while i < self.geno_names.len() && j < d.data.len() {
            if i != d.data[j] as usize {
                if write_index != i {
                    // Move the elements to their new positions
                    self.matrix_bit.swap(write_index, i);
                    self.geno_names.swap(write_index, i);
                    self.bim_entries.swap(write_index, i);
                    self.window_number.swap(write_index, i);
                }
                write_index += 1;
                if i < d.data[j] as usize{
                    i += 1;
                } else {
                    j += 1;
                }
            } else {
                i += 1;
                j += 1;
            }
        }
        if i != self.geno_names.len() {
            for x in i..self.geno_names.len() {
                if write_index != x {
                    // Move the elements to their new positions
                    self.matrix_bit.swap(write_index, x);
                    self.geno_names.swap(write_index, x);
                    self.bim_entries.swap(write_index, x);
                    self.window_number.swap(write_index, x);
                }
                write_index += 1;
            }
        }

        println!("Write index: {:?}", write_index);
        self.matrix_bit.truncate(write_index);
        self.geno_names.truncate(write_index);
        self.bim_entries.truncate(write_index);
        self.window_number.truncate(write_index);
    }

    //----------------------------------------------------------------------------------
    /// Write "empty" fam with no phenotypes
    ///
    /// Contains the names of the samples in the same order as plink bed file
    /// Stays the same for all runs
    pub fn write_fam(&self, number: usize, out_prefix: &str, _feature: Feature, len: usize) {
        let mut output = [out_prefix, &number.to_string(), "fam"].join(".");
        if len == 1 {
            output = [out_prefix, "fam"].join(".");
        }
        let f = File::create(output).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        if self.fam_entries.is_empty() {
            for x in self.sample_names.iter() {
                writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}", x, x, 0, 0, 0, 0)
                    .expect("Can not write file");
            }
        } else {
            for x in self.fam_entries.iter() {
                writeln!(f, "{}", x).expect("Can not write file");
            }
        }
    }

    /// Write bed file in SNP-major mode
    pub fn write_bed(&self, number: usize, out_prefix: &str, _feature: Feature, len: usize) {
        //hexdump -C test.bin
        // xxd -b file
        // xxd file
        // SNP: 00000001 , 0
        // IND: 00000000, 1

        let mut buff: Vec<u8> = vec![108, 27, 1];
        // Make SNP Vector

        for sel in self.matrix_bit.iter() {
            buff.extend(sel.as_raw_slice());
        }

        let mut output = [out_prefix, &number.to_string(), "bed"].join(".");
        if len == 1 {
            output = [out_prefix, "bed"].join(".");
        }
        let mut file = File::create(output).expect("Not able to write ");
        file.write_all(&buff).expect("Not able to write ")
    }

    pub fn write_chunks(&self, split: usize, output_prefix: &str, feature: Feature){
        let chunk_size = (self.matrix_bit.len() / split) + 1;
        let chunks = self.matrix_bit.chunks(chunk_size);

        let len = chunks.len();
        for (index, _y) in chunks.enumerate() {
            self.write_fam(index, output_prefix, feature, len);
            self.write_bed(index, output_prefix, feature, len);
            self.write_bim(index, output_prefix, &feature, len);
        }
    }


    /// Write bim file
    ///
    ///
    /// https://www.cog-genomics.org/plink/1.9/formats#bim
    /// Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
    /// Variant identifier
    /// Position in morgans or centimorgans (safe to use dummy value of '0')
    /// Base-pair coordinate (1-based; limited to 231-2)
    /// Allele 1 (corresponding to clear bits in .bed; usually minor)
    /// Allele 2 (corresponding to set bits in .bed; usually major)
    /// Representation here: [graph, ., 1, 0, A, T]
    pub fn write_bim(&self, number: usize, out_prefix: &str, feature: &Feature, len: usize) {
        let mut output = [out_prefix, &number.to_string(), "bim"].join(".");
        if len == 1 {
            output = [out_prefix, "bim"].join(".");
        }
        let file = File::create(output).expect("Unable to create file");
        let mut f = BufWriter::new(file);

        if self.bim_entries.is_empty() {
            for x in self.geno_names.iter() {
                writeln!(f, "graph\t.\t{}\t{}\tA\tT", 0, self.feature.to_string_u64(*x))
                    .expect("Can not write file");
            }

        } else {
            for x in self.bim_entries.iter() {
                writeln!(f, "{}", x).expect("Can not write file");
            }
        }
    }

    pub fn write_bimbam(
        &self,
        number: usize,
        out_prefix: &str,
        len: usize,
        val: &Vec<f32>,
    ) {
        let mut output = [out_prefix, &number.to_string(), "bimbam"].join(".");
        if len == 1 {
            output = [out_prefix, "bimbam"].join(".");
        }
        let file = File::create(output).expect("Unable to create file");
        let mut f = BufWriter::new(file);

        if self.bim_entries.is_empty() {
            for ((i, x), thresh) in self.geno_names.iter().enumerate().zip(val.iter()) {
                let p = normalize_vector(&self.matrix_u16[i], *thresh);
                let mut p2 = Vec::new();
                for x in self.sample_index_u16.iter() {
                    if x[0] == x[1] {
                        p2.push(p[x[0]]);
                    } else {
                        p2.push(average(&[p[x[0]], p[x[1]]]));
                    }
                }
                writeln!(
                    f,
                    "{}, A, T, {}",
                    self.feature.to_string_u64(*x),
                    p2.iter()
                        .map(|n| n.to_string())
                        .collect::<Vec<String>>()
                        .join(",  ")
                )
                .expect("Can not write file");
            }
        }
    }

    /// Write "empty" pheno file
    pub fn write_phenotype_bimbam(&self, number: usize, out_prefix: &str, len: usize) {
        let mut output = [out_prefix, &number.to_string(), "pheno"].join(".");
        if len == 1 {
            output = [out_prefix, "pheno"].join(".");
        }
        let f = File::create(output).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        if self.fam_entries.is_empty() {
            for x in self.sample_names.iter() {
                writeln!(f, "{}", x).expect("Can not write file");
            }
        }
    }
}

/// Normalize the vector for bimbam
fn normalize_vector(vector: &Vec<u16>, value: f32) -> Vec<f64> {
    let mut a = Vec::new();
    for element in vector.iter() {
        let result = (*element as f64 / value as f64).min(1.0);
        if result >= 1.0 {
            a.push(2.0)
        } else {
            a.push(result * 2.0);
        }
    }

    a
}

fn average(l: &[f64; 2]) -> f64 {
    (l[0] + l[1]) / 2.0
}
