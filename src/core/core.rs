use crate::core::helper::{is_all_ones, is_all_zeros, merge_u32_to_u64, to_string1, Feature};

use crate::r#mod::input_data::FileData;
use bitvec::prelude::*;
use gfa_reader::NCGfa;
use hashbrown::HashSet;
use std::fmt::Debug;
use std::fs::{metadata, File};
use std::io::{BufWriter, Write};
use std::ptr::write;

#[derive(Debug, Clone, Eq, PartialEq)]
/// Core data structure
///
/// SNP-major representation of the genotypes (also in memory) - Reading data takes longer here
/// Can represent two haplotypes (diploid) or one haplotype (haploid)
pub struct MatrixWrapper {
    pub shape: (usize, usize),
    pub matrix_core: Vec<Vec<u32>>,
    pub matrix_bin: Vec<BitVec<u8, Lsb0>>, //Vec<Vec<bool>>

    // Check if node, edges, dirnode, or alignment
    pub feature: Feature,

    pub geno_names: Vec<u64>,      // Name of all
    pub sample_names: Vec<String>, // Sample names (same order as in the matrix)

    // Plink specific stuff
    pub fam_entries: Vec<String>,
    pub bim_entries: Vec<String>,
}

impl MatrixWrapper {
    /// Dummy initialization
    pub fn new() -> Self {
        Self {
            shape: (0, 0),
            matrix_core: Vec::new(),
            matrix_bin: Vec::new(),
            feature: Feature::Node,

            // SNP
            geno_names: Vec::new(),
            bim_entries: Vec::new(),

            // Fam
            sample_names: Vec::new(),
            fam_entries: Vec::new(),
        }
    }

    /// Initialize the SNP index
    pub fn make_index(&mut self, data: &NCGfa<()>, t: Feature) {
        //let mut bb = HashMap::with_hasher(BuildHasherDefault::<NoHashHasher<u32>>::default());
        let mut geno_names = Vec::new();
        match t {
            Feature::Node => {
                for (i, x) in data.nodes.iter().enumerate() {
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

                    for (i, x) in edd2.iter().enumerate() {
                        geno_names.push(*x);
                    }
                }
            }
            Feature::Edge => {
                if data.edges.is_some() {
                    let value = data.edges.as_ref().unwrap();
                    for (i, x) in value.iter().enumerate() {
                        let (v1, v2, v3, v4) = (x.from, x.from_dir, x.to, x.to_dir);
                        let u1 = v1 * 2 + v2 as u32;
                        let u2 = v3 * 2 + v4 as u32;
                        let uu = merge_u32_to_u64(u1, u2);
                        geno_names.push(uu);
                    }
                }
            }
        }
        self.geno_names = geno_names;
    }

    pub fn remove_non_info(&mut self) {
        let mut write_index = 0;

        for read_index in 0..self.geno_names.len() {
            if !(is_all_ones(&self.matrix_bin[read_index])
                || is_all_zeros(&self.matrix_bin[read_index]))
            {
                // Retain elements that satisfy the condition
                if write_index != read_index {
                    // Move the elements to their new positions
                    self.matrix_bin.swap(write_index, read_index);
                    self.geno_names.swap(write_index, read_index);
                }

                write_index += 1;
            }
        }
        self.matrix_bin.truncate(write_index);
        self.geno_names.truncate(write_index);
    }

    pub fn remove_paths(&mut self, to_be_removed: &Vec<String>) {
        let mut remove_index: Vec<usize> = Vec::new();
        let remo: HashSet<String> = to_be_removed.iter().cloned().collect();
        for (i, x) in self.sample_names.iter().enumerate() {
            if remo.contains(x) {
                remove_index.push(i);
            }
        }


        self.sample_names.retain(|x| !remo.contains(x));
        if !self.fam_entries.is_empty() {
            self.fam_entries.retain(|x| {
                !remo.contains(&x.split_whitespace().collect::<Vec<&str>>()[0].to_string())
            });
        }
        assert_eq!(self.sample_names.len(), self.fam_entries.len());
        let mut i = 0;
        while i < self.matrix_bin.len() {
            let mut a = &mut self.matrix_bin[i];
            let mut f = 0;
            for x in remove_index.iter() {
                a.remove(*x * 2 - f * 2);
                a.remove(*x * 2 - f * 2 + 1);
                f += 1
            }
            i += 1;
        }
    }

    pub fn remove_feature(&mut self, d: &FileData) {
        let mut write_index = 0;
        let mut i = 0;
        let mut j = 0;
        while i < self.geno_names.len() && j < d.data.len() {
            if self.geno_names[i] != d.data[j] {
                if write_index != i {
                    // Move the elements to their new positions
                    self.matrix_bin.swap(write_index, i);
                    self.geno_names.swap(write_index, i);
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
        self.matrix_bin.truncate(write_index);
        self.geno_names.truncate(write_index);
    }

    /// Write "empty" fam with no phenotypes
    ///
    /// Contains the names of the samples in the same order as plink bed file
    pub fn write_fam(&self, number: usize, out_prefix: &str, feature: Feature, len: usize) {
        let mut output = [out_prefix, &feature.to_string1(), &number.to_string(), "fam"].join(".");
        if len == 1 {
            output = [out_prefix, &feature.to_string1(), "fam"].join(".");
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
    pub fn write_bed(&self, number: usize, out_prefix: &str, feature: Feature, len: usize) {
        //hexdump -C test.bin
        // xxd -b file
        // xxd file
        // SNP: 00000001 , 0
        // IND: 00000000, 1

        let mut buff: Vec<u8> = vec![108, 27, 1];
        // Make SNP Vector

        for sel in self.matrix_bin.iter() {
            buff.extend(sel.as_raw_slice());
        }

        let mut output = [out_prefix, &feature.to_string1(), &number.to_string(), "bed"].join(".");
        if len == 1 {
            output = [out_prefix, &feature.to_string1(), "bed"].join(".");
        }
        let mut file = File::create(output).expect("Not able to write ");
        file.write_all(&buff).expect("Not able to write ")
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
        let mut output = [
            out_prefix,
            &feature.to_string1(),
            &number.to_string(),
            "bim",
        ]
        .join(".");
        if len == 1 {
            output = [out_prefix, &feature.to_string1(), "bim"].join(".");
        }
        let f = File::create(output).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        if self.bim_entries.is_empty() {
            for x in self.geno_names.iter() {
                writeln!(f, "graph\t.\t{}\t{}\tA\tT", 0, to_string1(*x, feature))
                    .expect("Can not write file");
            }
        } else {
            for x in self.bim_entries.iter() {
                writeln!(f, "{}", x).expect("Can not write file");
            }
        }
    }
}
