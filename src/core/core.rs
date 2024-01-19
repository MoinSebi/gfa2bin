use crate::core::helper::{is_all_ones, is_all_zeros, Feature, GenoName};

use bitvec::prelude::*;
use gfa_reader::NCGfa;
use nohash_hasher::NoHashHasher;
use std::collections::HashMap;
use std::fmt::Debug;
use std::fs::File;
use std::hash::BuildHasherDefault;
use std::io::{BufWriter, Write};

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

    pub geno_names: Vec<GenoName>, // Name of all
    pub geno_map: HashMap<GenoName, usize, BuildHasherDefault<NoHashHasher<u32>>>, // Mapping from "name" to index in the matrix
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

            geno_names: Vec::new(),
            geno_map: HashMap::default(),
            sample_names: Vec::new(),

            fam_entries: Vec::new(),
            bim_entries: Vec::new(),
        }
    }

    /// Initialize the SNP index
    pub fn make_index(&mut self, data: &NCGfa<()>, t: Feature) {
        match t {
            Feature::Node => {
                for (i, x) in data.nodes.iter().enumerate() {
                    self.geno_map.insert(GenoName { name: x.id as u64 }, i);
                    self.geno_names.push(GenoName { name: x.id as u64 });
                }
            }
            Feature::DirNode => {
                if data.edges.is_some() {
                    let value = data.edges.as_ref().unwrap();
                    for (i, x) in value.iter().enumerate() {
                        self.geno_map.insert(
                            GenoName {
                                name: x.from as u64 * 2 + x.from_dir as u64,
                            },
                            i,
                        );
                        self.geno_names.push(GenoName {
                            name: x.from as u64 * 2 + x.from_dir as u64,
                        });
                    }
                }
            }
            Feature::Edge => {
                if data.edges.is_some() {
                    let value = data.edges.as_ref().unwrap();
                    for (i, x) in value.iter().enumerate() {
                        self.geno_map.insert(
                            GenoName {
                                name: x.from as u64 * 2 + x.from_dir as u64,
                            },
                            i,
                        );
                        self.geno_names.push(GenoName {
                            name: x.from as u64 * 2 + x.from_dir as u64,
                        });
                    }
                }
            }
        }
    }

    pub fn remove_non_info(&mut self) {
        let mut write_index = 0;
        let mut remove_index = Vec::new();

        for read_index in 0..self.geno_names.len() {
            if !is_all_ones(&self.matrix_bin[read_index])
                || !is_all_zeros(&self.matrix_bin[read_index])
            {
                println!("dsajkdja {:?}", self.matrix_bin[read_index]);
                remove_index.push(read_index);
                // Retain elements that satisfy the condition
                if write_index != read_index {
                    // Move the elements to their new positions
                    self.matrix_bin.swap(write_index, read_index);
                    self.sample_names.swap(write_index, read_index);
                }

                write_index += 1;
            } else {
                println!("dsajkdasdadadsadja {:?}", self.matrix_bin[read_index]);
            }
        }
        println!("{:?}", write_index);
        self.matrix_bin.truncate(write_index);
        self.geno_names.truncate(write_index);

        let mut k = self
            .geno_map
            .iter()
            .map(|(k, v)| (*v, k.clone()))
            .collect::<Vec<(usize, GenoName)>>();
        k.sort_by(|a, b| a.0.cmp(&b.0));
        let mut gg = 0;
        for x in remove_index {
            for y in gg..k.len() {
                if k[y].0 == x {
                    self.geno_map.remove(&k[y].1);
                    gg = y;
                    break;
                }
            }
            self.geno_map.remove(&k[x].1);
        }
    }

    /// Write "empty" fam with no phenotypes
    ///
    /// Contains the names of the samples in the same order as plink bed file
    pub fn write_fam(&self, number: usize, out_prefix: &str, feature: &str, len: usize) {
        let mut output = [out_prefix, feature, &number.to_string(), "fam"].join(".");
        if len == 1 {
            output = [out_prefix, feature, "fam"].join(".");
        }
        let f = File::create(output).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        for x in self.sample_names.iter() {
            writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}", x, x, 0, 0, 0, 0).expect("Can not write file");
        }
    }

    /// Write bed file in SNP-major mode
    pub fn write_bed(&self, number: usize, out_prefix: &str, feature: &str, len: usize) {
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

        let mut output = [out_prefix, feature, &number.to_string(), "bed"].join(".");
        if len == 1 {
            output = [out_prefix, feature, "bed"].join(".");
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
        let mut output = [out_prefix, &feature.to_string(), &number.to_string(), "bim"].join(".");
        if len == 1 {
            output = [out_prefix, &feature.to_string(), "bim"].join(".");
        }
        let f = File::create(output).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        for x in self.geno_names.iter() {
            writeln!(f, "graph\t.\t{}\t{}\tA\tT", 0, x.to_string(feature))
                .expect("Can not write file");
        }
    }

    pub fn filter_shared2(&mut self) {
        let mut write_index = 0;
        for read_index in 0..self.sample_names.len() {
            if is_all_ones(&self.matrix_bin[read_index])
                || is_all_zeros(&self.matrix_bin[read_index])
            {
                // Retain elements that satisfy the condition
                if write_index != read_index {
                    // Move the elements to their new positions
                    self.matrix_bin.swap(write_index, read_index);
                    self.sample_names.swap(write_index, read_index);
                }

                write_index += 1;
            }
        }
        let a = &self.sample_names[0..write_index];

        for x in a {
            self.geno_map.remove(&GenoName {
                name: x.parse::<u64>().unwrap(),
            });
        }
        // Truncate both vectors to their new length
        self.matrix_bin.truncate(write_index);
        self.sample_names.truncate(write_index);
    }
}
