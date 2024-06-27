use crate::core::helper::{
    is_all_ones, is_all_zeros, merge_u32_to_u64, Feature,
};

use crate::r#mod::input_data::FileData;
use crate::r#mod::mod_main::remove_feature;
use bitvec::prelude::*;
use gfa_reader::Gfa;
use hashbrown::HashSet;
use log::info;


use std::fmt::Debug;
use std::fs::File;
use std::io::{BufWriter, Write};

#[derive(Debug, Clone, PartialEq)]
/// Core data structure
///
/// SNP-major representation of the genotypes (also in memory) - Reading data takes longer here
/// Can represent two haplotypes (diploid) or one haplotype (haploid)
pub struct MatrixWrapper {
    pub shape: (usize, usize),             //Update
    pub matrix_u16: Vec<Vec<u16>>,         // Raw values
    pub matrix_f32: Vec<Vec<f32>>,         // Normalized values
    pub matrix_bit: Vec<BitVec<u8, Lsb0>>, //Vec<Vec<bool>>

    // Check if node, edges, dirnode, or alignment
    pub feature: Feature,          // Feature - DIRNODE, NODE, EDGE or mode
    pub geno_names: Vec<u64>,      // Name of all - "SNP" names
    pub window_number: Vec<u32>,   // Number of window
    pub window_size: usize,        // Size of windows
    pub sample_names: Vec<String>, // Sample names (same order as in the matrix)
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
            matrix_f32: Vec::new(),
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

    /// Initialize the "SNP" index
    /// Index is used to
    pub fn make_index(&mut self, data: &Gfa<u32, (), ()>, t: Feature) {
        //let mut bb = HashMap::with_hasher(BuildHasherDefault::<NoHashHasher<u32>>::default());
        let mut geno_names = Vec::new();
        match t {
            Feature::Node => {
                for (_i, x) in data.segments.iter().enumerate() {
                    geno_names.push(x.id as u64);
                    //bb.insert(x.id as u64, i);
                }
            }
            Feature::DirNode => {
                if !data.links.is_empty() {
                    let value = &data.links;
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
                if !data.links.is_empty() {
                    let value = &data.links;
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

    /// Create a presence/absence matrix based on a threshold
    pub fn matrix2bin(&mut self, relative: &Vec<f32>) {
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

    fn indices_where<F, T>(&self, condition: F, data: Vec<T>) -> Vec<usize>
    where
        F: Fn(&T) -> bool,
    {
        data.iter()
            .enumerate()
            .filter_map(|(index, item)| if condition(item) { Some(index) } else { None })
            .collect()
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

    pub fn remove_samples2(&mut self, to_be_removed: &Vec<usize>) {
        let mut i = 0;
        while i < self.matrix_bit.len() {
            let a = &mut self.matrix_bit[i];
            let mut f = 0;
            for x in to_be_removed.iter() {
                a.remove(*x * 2 - f * 2 + 1);
                a.remove(*x * 2 - f * 2);
                f += 1;
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

    pub fn remove_test1(&self, file_data: FileData) -> Vec<usize> {
        if file_data.feature == Feature::MWindow || file_data.feature == Feature::Block {
            let bb = self.geno_names.iter().zip(self.window_number.iter());
            let cc = file_data.data.iter().zip(file_data.window.iter());

            
            remove_feature(bb, cc)
        } else {
            let bb = self.geno_names.iter();
            let cc = file_data.data.iter();

            
            remove_feature(bb, cc)
        }
    }

    /// Remove entries by index from the matrix
    ///
    /// Removed in:
    /// - genome names
    /// - matrix bit
    /// - bim entries
    ///
    pub fn remove_by_index(&mut self, index: Vec<usize>) {
        let mut rmi = self.matrix_bit.len() - 1;
        for x in index.iter().rev() {
            self.matrix_bit.swap(rmi, *x);
            self.geno_names.swap(rmi, *x);
            self.bim_entries.swap(rmi, *x);
            rmi -= 1;
        }
        if !self.window_number.is_empty() {
            for x in index.iter().rev() {
                self.window_number.swap(*x, rmi);
                rmi -= 1;
            }
        }

        self.matrix_bit.truncate(rmi);
        self.geno_names.truncate(rmi);
        self.window_number.truncate(rmi);
        self.bim_entries.truncate(rmi);
    }

    //----------------------------------------------------------------------------------
    /// Write wrapper
    pub fn write_wrapper(
        &self,
        b: bool,
        split: usize,
        output_prefix: &str,
        thresh: Vec<f32>,
        feature_enum: Feature,
    ) {
        if b {
            info!("Writing the bimbam");
            let chunk_size = (self.matrix_u16.len() / split) + 1;
            let chunks = self.matrix_u16.chunks(chunk_size);
            let len = chunks.len();

            for (index, _y) in chunks.enumerate() {
                self.write_bimbam(index, output_prefix, len, &thresh);
                self.write_phenotype_bimbam(index, output_prefix, len);
            }
        } else {
            // Output
            info!("Writing the output");
            let chunk_size = (self.matrix_bit.len() / split) + 1;
            let chunks = self.matrix_bit.chunks(chunk_size);

            let len = chunks.len();
            for (index, _y) in chunks.enumerate() {
                //write_bed2(y, output_prefix, feature, index, len);
                self.write_fam(index, output_prefix, feature_enum, len);
                self.write_bed(index, output_prefix, feature_enum, len);
                self.write_bim(index, output_prefix, &feature_enum, len);
            }
        }
    }

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
                writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}", x, x, 0, 0, 0, 1)
                    .expect("Can not write file");
            }
        } else {
            for x in self.fam_entries.iter() {
                writeln!(f, "{}", x).expect("Can not write file");
            }
        }
    }

    /// Write bed file in SNP-major mode
    ///
    /// https://zzz.bwh.harvard.edu/plink/binary.shtml
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

    /// Write chunks (splits)
    ///
    /// Split the data into chunks and write them
    /// - bed
    /// - bim
    /// - fam
    pub fn write_chunks(&self, split: usize, output_prefix: &str, feature: Feature) {
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
    pub fn write_bim(&self, number: usize, out_prefix: &str, _feature: &Feature, len: usize) {
        let mut output = [out_prefix, &number.to_string(), "bim"].join(".");
        if len == 1 {
            output = [out_prefix, "bim"].join(".");
        }
        let file = File::create(output).expect("Unable to create file");
        let mut f = BufWriter::new(file);

        if self.bim_entries.is_empty() {
            for x in self.geno_names.iter() {
                writeln!(
                    f,
                    "graph\t.\t{}\t{}\tA\tT",
                    0,
                    self.feature.to_string_u64(*x)
                )
                .expect("Can not write file");
            }
        } else {
            for x in self.bim_entries.iter() {
                writeln!(f, "{}", x).expect("Can not write file");
            }
        }
    }

    /// Write a bimbam file
    ///
    /// Based on real values (no presence/absence) and a threshold
    /// Default genotype is A and T
    pub fn write_bimbam(&self, number: usize, out_prefix: &str, len: usize, val: &Vec<f32>) {
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

    /// Write "empty" pheno file for bimbam
    ///
    /// Format - simple list of sample names
    /// - Sample name
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
