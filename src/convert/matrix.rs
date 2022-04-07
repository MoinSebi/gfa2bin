
use std::fmt::{Debug};
use std::fs::File;
use std::io::{Write, BufWriter};
use log::info;
use crate::matrix_edge;


/// Core structure
/// 2D Matrix with
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub struct Matrix {
    pub shape: (u32, u32),
    pub matrix_core: Vec<Vec<u32>>,
}


impl Matrix {
    /// Constructor
    pub fn new() -> Self {
        let shape: (u32, u32) = (0,0);
        let matrix: Vec<Vec<u32>> = Vec::new();
        Self {
            shape: shape,
            matrix_core: matrix,
        }
    }

    /// Make binary
    /// @param:
    ///
    /// Input: Threshold
    /// Basic: number = threshhold
    /// For presence/absence -> number = 1
    pub fn copy(&self, number: u32) -> Vec<bitvec::vec::BitVec>{
        let mut new_matrix: Vec<bitvec::vec::BitVec> = Vec::new();
        for x in self.matrix_core.iter(){
            let mut new_vec: bitvec::vec::BitVec = bitvec::vec::BitVec::new();
            for y in x.iter(){
                new_vec.push(y.clone() >= number);

            }
            new_matrix.push(new_vec);
        }
        new_matrix
    }

    #[allow(dead_code)]
    /// MinMaxScaler
    pub fn min_max(&self) -> Vec<Vec<f32>>{
        let mut new_matrix:  Vec<Vec<f32>> = Vec::new();
        for x in 0..self.matrix_core[0].len(){
            let mut max = 0;
            for y in 0..self.matrix_core.len(){
                if self.matrix_core[y][x] > max{
                    max = self.matrix_core[y][x]
                }
            }
            let mut vec_new: Vec<f32> = Vec::new();
            for y in 0..self.matrix_core.len(){
                vec_new.push(self.matrix_core[y][x] as f32/max as f32) //this is wrong
            }
            new_matrix.push(vec_new);
        }

        new_matrix
    }


    pub fn write_bimbam(&self, out_prefix: &str, t: &str){
        info!("Writing bimbam file");
        let f = File::create([out_prefix, t, "bed"].join(".")).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        let min_max_matrix = self.min_max();
        for (index, x) in min_max_matrix.iter().enumerate(){
            let j: Vec<String> = x.iter().map(|i| format!("{}", i)).collect();
            write!(f, "{}\t{}\t{}\t", index, "A", "T").expect("Not able to write ");
            write!(f, "{}\n", j.join("\t")).expect("Not able to write");
        }
    }




}




















