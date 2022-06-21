
use std::fmt::{Debug};
use std::fs::File;
use std::io::{Write, BufWriter};
use log::info;

/// Core structure
/// 2D Matrix with
#[derive(Debug, Clone, Eq, PartialEq, Hash)]
pub struct Matrix {
    pub shape: (u32, u32),
    pub matrix_core: Vec<Vec<u32>>,
}


impl Matrix {
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

    #[allow(dead_code)]
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




















