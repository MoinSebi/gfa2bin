use crate::core::helper::{merge_u32_to_u64, Feature};


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
    // Matrix
    pub shape: (usize, usize),             // Update
    pub matrix_u16: Vec<Vec<u16>>,         // Raw values
    pub matrix_f32: Vec<Vec<f32>>,         // Normalized values
    pub matrix_bit: Vec<BitVec<u8, Lsb0>>, // BitVec Vec<Vec<bool>>

    // SNP
    pub feature: Feature, // Feature - DIRNODE, NODE, EDGE, subpath, block, window
    pub window_number: Vec<u32>, // Number of window
    pub window_size: usize, // Size of windows
    pub geno_names: Vec<u64>, // Name of all - "SNP" names
    pub bim_entries: Vec<String>, // Bim entries

    // Fam - Samples
    pub sample_names: Vec<String>, // Sample names (same order as in the matrix)
    pub fam_entries: Vec<String>,  // Fam entries
    pub sample_index_u16: Vec<[usize; 2]>, // Sample index [11, 12] =
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
    ///
    /// The SNP index is used for fast insertion of data
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

    pub fn matrix2bin<T>(
        input_data: &Vec<Vec<T>>,
        relative: &Vec<f32>,
        sample_index: &Vec<[usize; 2]>,
    ) -> Vec<BitVec<u8>>
    where
        T: PartialOrd + Copy + Into<f64>,
    {
        let mut matrix_bin = Vec::new();
        for (val, re) in input_data.iter().zip(relative.iter()) {
            let mut biit = BitVec::new();
            for aa in sample_index.iter() {
                if aa[0] == aa[1] {
                    let a = val[aa[0]].into();
                    if a >= *re as f64 {
                        biit.push(true);
                        biit.push(true);
                    } else {
                        biit.push(false);
                        biit.push(false);
                    }
                } else {
                    let a = val[aa[0]].into();
                    let b = val[aa[1]].into();

                    let a1 = a >= *re as f64;
                    let b1 = b >= *re as f64;

                    if a1 && b1 {
                        biit.push(true);
                        biit.push(true);
                    } else if a1 != b1 {
                        biit.push(false);
                        biit.push(true);
                    } else {
                        biit.push(false);
                        biit.push(false);
                    }
                }
            }
            matrix_bin.push(biit)
        }
        matrix_bin
    }

    //----------------------------------------------------------------------------------
    /// Write wrapper
    pub fn write_wrapper(
        &mut self,
        bimbam: bool,
        split: usize,
        output_prefix: &str,
        thresh: Vec<f32>,
        feature_enum: Feature,
        pheno: f64,
        remove_non_info: bool,
    ) {
        if bimbam {
            info!("Writing the bimbam");

            if self.matrix_f32.is_empty() {
                let chunk_size = (self.matrix_u16.len() / split) + 1;
                let chunks = self.matrix_u16.chunks(chunk_size);
                let len = chunks.len();
                for (index, y) in chunks.enumerate() {
                    self.write_bimbam(index, output_prefix, len, &thresh, y);
                    self.write_phenotype_bimbam(index, output_prefix, len, pheno);
                }
            } else {
                let chunk_size = (self.matrix_f32.len() / split) + 1;
                let chunks = self.matrix_f32.chunks(chunk_size);
                let len = chunks.len();
                for (index, y) in chunks.enumerate() {
                    self.write_bimbam(index, output_prefix, len, &thresh, y);
                    self.write_phenotype_bimbam(index, output_prefix, len, pheno);
                }
            }

            // if plink bed
        } else {
            if self.matrix_bit.is_empty() {
                if !self.matrix_f32.is_empty() {
                    self.matrix_bit = MatrixWrapper::matrix2bin(
                        &self.matrix_f32,
                        &thresh,
                        &self.sample_index_u16,
                    );
                } else {
                    self.matrix_bit = MatrixWrapper::matrix2bin(
                        &self.matrix_u16,
                        &thresh,
                        &self.sample_index_u16,
                    );
                }
            }
            info!(
                "Matrix [SNPs X Samples] (before remove): {}, {}",
                self.matrix_bit.len(),
                self.matrix_bit[0].len()
            );
            if remove_non_info {
                //self.remove_non_info();
                info!(
                    "Matrix [SNPs X Samples] (after remove): {}, {}",
                    self.matrix_bit.len(),
                    self.matrix_bit[0].len()
                );
            }

            // Output
            info!("Writing the plink bed/bim/fam");
            self.write_chunks(split, output_prefix, feature_enum, pheno);
        }
    }

    /// Write chunks (splits)
    ///
    /// Split the data into chunks and write them
    /// - bed
    /// - bim
    /// - fam
    pub fn write_chunks(&self, split: usize, output_prefix: &str, feature: Feature, pheno: f64) {
        let chunk_size = (self.matrix_bit.len() / split) + 1;
        let chunks = self.matrix_bit.chunks(chunk_size);

        let len = chunks.len();
        for (index, _y) in chunks.enumerate() {
            self.write_fam(index, output_prefix, feature, len, pheno);
            self.write_bed(index, output_prefix, feature, len);
            self.write_bim(index, output_prefix, &feature, len);
        }
    }

    /// Write "empty" fam with no phenotypes
    ///
    /// Contains the names of the samples in the same order as plink bed file
    /// Stays the same for all runs
    pub fn write_fam(
        &self,
        number: usize,
        out_prefix: &str,
        _feature: Feature,
        len: usize,
        mut pheno: f64,
    ) {
        let mut output = [out_prefix, &number.to_string(), "fam"].join(".");
        if len == 1 {
            output = [out_prefix, "fam"].join(".");
        }
        if pheno == f64::MAX {
            pheno = -9.0
        }

        let f = File::create(output).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        if self.fam_entries.is_empty() {
            for x in self.sample_names.iter() {
                writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}", x, x, 0, 0, 0, pheno)
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
                writeln!(f, "graph\t.\t{}\t{}\tA\tT", 0, feature.to_string_u64(*x))
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
    pub fn write_bimbam<T>(
        &self,
        number: usize,
        out_prefix: &str,
        len: usize,
        val: &Vec<f32>,
        vv: &[Vec<T>],
    ) where
        T: Into<f64> + Copy,
    {
        let mut output = [out_prefix, &number.to_string(), "bimbam"].join(".");
        if len == 1 {
            output = [out_prefix, "bimbam"].join(".");
        }
        let file = File::create(output).expect("Unable to create file");
        let mut f = BufWriter::new(file);

        if self.bim_entries.is_empty() {
            for ((i, x1), thresh) in self.geno_names.iter().enumerate().zip(val.iter()) {
                let p = normalize_vector(&vv[i], *thresh as f64);
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
                    self.feature.to_string_u64(*x1),
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
    pub fn write_phenotype_bimbam(&self, number: usize, out_prefix: &str, len: usize, pheno: f64) {
        let mut output = [out_prefix, &number.to_string(), "pheno"].join(".");
        if len == 1 {
            output = [out_prefix, "pheno"].join(".");
        }

        let mut pheno_string = pheno.to_string();
        if pheno == f64::MAX {
            pheno_string = "NA".to_string();
        }

        let f = File::create(output).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        if self.fam_entries.is_empty() {
            for _x in self.sample_names.iter() {
                writeln!(f, "{}", pheno).expect("Can not write file");
            }
        } else {
            panic!("Text")
        }
    }
}

/// Normalize the vector
fn normalize_vector<T>(vector: &[T], value: f64) -> Vec<f64>
where
    T: Into<f64> + Copy,
{
    let mut normalized = Vec::with_capacity(vector.len());

    for &element in vector {
        let value_f64 = element.into();
        let result = (value_f64 / value).min(1.0);
        let normalized_value = if result >= 1.0 { 2.0 } else { result * 2.0 };
        normalized.push(normalized_value);
    }

    normalized
}

fn average(l: &[f64; 2]) -> f64 {
    (l[0] + l[1]) / 2.0
}
