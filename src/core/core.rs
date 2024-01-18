use crate::core::helper::{Feature, GenoName};

use bitvec::prelude::*;
use gfa_reader::{NCGfa};
use std::collections::HashMap;
use std::fmt::{Debug};
use std::fs::{File};
use std::io::{BufWriter, Write};

#[derive(Debug, Clone, Eq, PartialEq)]
/// Core data structure
///
/// SNP-major representation of the genotypes
pub struct MatrixWrapper {
    pub shape: (usize, usize),
    pub matrix_core: Vec<Vec<u32>>,
    pub matrix_bin: Vec<BitVec<u8, Lsb0>>, //Vec<Vec<bool>>

    // Check if node, edges, dirnode, or alignment
    pub feature: Feature,

    pub geno_names: Vec<GenoName>,          // Name of all
    pub geno_map: HashMap<GenoName, usize>, // Mapping from "name" to index in the matrix
    pub sample_names: Vec<String>,          // Sample names (same order as in the matrix)

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
            geno_map: HashMap::new(),
            sample_names: Vec::new(),

            fam_entries: Vec::new(),
            bim_entries: Vec::new(),
        }
    }

    //--------------------------------------------------------------------------------

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

    // pub fn make_binary(& mut self, threshold: u32, haplotype: &Haplotype){
    //     let mut new_matrix =  vec![BitVec::repeat(false,haplotype.haplotype.len()*2); self.matrix_core.len()];
    //     for (i, x) in haplotype.haplotype2.iter().enumerate() {
    //         self.bin_names.push(x.0.clone());
    //         for (i2, y) in self.matrix_core.iter().enumerate() {
    //
    //             if y[x.1[0]] >= threshold {
    //                 new_matrix.get_mut(i2).unwrap().set(i*2 as usize, true);
    //             }
    //             if y[x.1[1]] >= threshold {
    //                 new_matrix.get_mut(i2).unwrap().set(i*2+1 as usize, true);
    //             }
    //         }
    //     }
    //     self.matrix_bin = new_matrix;
    // }

    //--------------------------------------------------------------------------------
    // Modification

    /// Write "empty" fam with no phenotypes
    /// Contains the names of the individuals in the same order as plink bed file
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

    pub fn write_bim(&self, number: usize, out_prefix: &str, feature: &Feature, len: usize) {
        let mut output = [out_prefix, &feature.to_string(), &number.to_string(), "bim"].join(".");
        if len == 1 {
            output = [out_prefix, &feature.to_string(), "bim"].join(".");
        }
        let f = File::create(output).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        for x in self.geno_names.iter() {
            writeln!(
                f,
                "graph\t.\t{}\t{}\tA\tT",
                0,
                x.to_string(feature)
            )
            .expect("Can not write file");
        }
    }

    // pub fn filter_shared2(&mut self){
    //     let mut write_index = 0;
    //     for read_index in 0..self.sample_names.len() {
    //         if is_all_ones(&self.matrix_bin[read_index])  || is_all_zeros(&self.matrix_bin[read_index]) {
    //             // Retain elements that satisfy the condition
    //             if write_index != read_index {
    //                 // Move the elements to their new positions
    //                 self.matrix_bin.swap(write_index, read_index);
    //                 self.sample_names.swap(write_index, read_index);
    //             }
    //
    //             write_index += 1;
    //         }
    //     }
    //     let a = &self.sample_names[0..write_index];
    //
    //     for x in a{
    //         self.geno_map.remove(&GenoName{name: x.parse::<u64>().unwrap()});
    //     }
    //     // Truncate both vectors to their new length
    //     self.matrix_bin.truncate(write_index);
    //     self.sample_names.truncate(write_index);
    //
    // }
}
