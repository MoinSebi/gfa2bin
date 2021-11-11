
use std::fmt::{Debug};
use std::fs::File;
use std::io::{Write, BufWriter};
use crate::helper::{transpose, trans2};
use std::collections::HashSet;


/// Core structure
/// 2D Matrix with
#[derive(Debug, Clone)]
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
    ///
    /// Input: Threshold
    /// Basic: number = threshhold
    /// For presence/absence -> number = 1
    pub fn copy(&self, number: u32) -> Vec<Vec<bool>>{
        let mut new_matrix: Vec<Vec<bool>> = Vec::new();
        for x in self.matrix_core.iter(){
            let mut new_vec: Vec<bool> = Vec::new();
            for y in x.iter(){
                new_vec.push(y.clone() >= number);

            }
            new_matrix.push(new_vec);
        }
        new_matrix
    }


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


    /// Write bimbam file
    pub fn write_bimbam(&self, out_prefix: &str, t: &str){
        println!("Writing bimbam file");
        let f = File::create([out_prefix, t, "bed"].join(".")).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        let min_max_matrix = self.min_max();
        for (index, x) in min_max_matrix.iter().enumerate(){
            let j: Vec<String> = x.iter().map(|i| format!("{}", i)).collect();
            write!(f, "{}\t{}\t{}\t", index, "A", "T").expect("Not able to write ");
            write!(f, "{}\n", j.join("\t")).expect("Not able to write");
        }
    }

    pub fn filter(&self) -> Vec<usize>{
        eprintln!("Filtering");
        println!("{} {}", self.matrix_core.len(), self.matrix_core[0].len());
        let k: Vec<Vec<u32>>= trans2(&self.matrix_core);
        let mut k2 = Vec::new();
        let mut count = 0;
        let mut kk: Vec<usize> = Vec::new();

        for (i, x) in k.iter().enumerate(){
            let mut sum = 0;
            for y in x.iter() {
                if y != &0 {
                    sum += 1;
                }
            }
            if sum as usize != x.len(){
                k2.push(x.clone());
            } else {
                println!("{} {}", sum, x.len());
                kk.push(i);
                count += 1;
            }


        }
        eprintln!("{}", kk.len());
        let k3 = transpose(&k2);

        eprintln!("Before {}  After {}", self.matrix_core[0].len(), k3[0].len());
        return kk;
    }

    pub fn reduce_comb(&self){
        let mut hs: HashSet<_> = HashSet::new();
        let k: Vec<Vec<u32>>= transpose(&self.matrix_core);
        for x in k.iter() {
            hs.insert(x);
        }
        println!("Reduce {}", hs.len());
        println!("Reduce {:?}", hs);
    }



}




















