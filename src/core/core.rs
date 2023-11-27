

use std::collections::{HashMap};
use std::fmt::Debug;
use std::io::{Write, BufWriter};
use std::fs::File;
use bitvec::prelude::*;
use crate::graph::parser::Haplotype;
use crate::helper::{custom_retain_two_vectors, is_all_ones, is_all_zeros};


#[derive(Debug, Clone, Eq, PartialEq)]
/// Core data structure, which includes ever
pub struct MatrixWrapper {
    pub shape: (usize, usize),
    pub matrix_core: Vec<Vec<u32>>,
    pub core_names: HashMap<u32, String>,
    pub matrix_bin: Vec<BitVec<u8, Msb0>>, //Vec<Vec<bool>>
    pub transposed: bool,
    pub names: Vec<String>,
    pub bin_names: Vec<String>,
    pub haplo: Vec<(String, Vec<usize>)>,

}



impl MatrixWrapper {

    /// Dummy initialization
    pub fn new() -> Self {
        let col: HashMap<u32, String> = HashMap::new();
        let bv: Vec<BitVec<u8, Msb0>> = Vec::new();
        let matrix = Vec::new();
        Self {
            shape: (0,0),
            matrix_core: matrix,
            core_names: col,
            matrix_bin: bv,
            transposed: false,
            names: Vec::new(),
            haplo: Vec::new(),
            bin_names: Vec::new(),
        }
    }


    //--------------------------------------------------------------------------------


    pub fn make_binary2(& mut self, number: u32, haplotype: &Haplotype){
        let mut new_matrix =  vec![BitVec::repeat(false,haplotype.haplotype.len()*2); self.matrix_core.len()];
        for (i, x) in haplotype.haplotype2.iter().enumerate() {
            self.bin_names.push(x.0.clone());
            for (i2, y) in self.matrix_core.iter().enumerate() {

                if y[x.1[0]] >= number {
                    new_matrix.get_mut(i2).unwrap().set(i*2 as usize, true);
                }
                if y[x.1[1]] >= number {
                    new_matrix.get_mut(i2).unwrap().set(i*2+1 as usize, true);
                }
            }
        }
        self.matrix_bin = new_matrix;
    }



    //--------------------------------------------------------------------------------
    // Modification













    /// Write "empty" fam with no phenotypes
    /// Contains the names of the individuals in the same order as plink bed file
    pub fn write_fam(&self, number: usize, out_prefix: &str, feature: &str, len: usize){
        let mut output = [out_prefix,  feature, &number.to_string(), "fam"].join(".");
        if len == 1{
            output = [out_prefix,  feature, "fam"].join(".");
        }
        let f = File::create(output).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        for x in self.bin_names.iter(){
            write!(f,"{}\t{}\t{}\t{}\t{}\t{}\n", x, x, 0, 0, 0, 0).expect("Can not write file");
        }
    }









    /// Filter out all entries which have
    /// - only zeros
    /// - only ones
    pub fn filter_shared1<T>(&mut self, vec1: &mut Vec<T>){
        custom_retain_two_vectors(&mut self.matrix_bin, vec1, |x| !is_all_zeros(x));
        custom_retain_two_vectors(&mut self.matrix_bin, vec1, |x| !is_all_ones(x));

    }



}


