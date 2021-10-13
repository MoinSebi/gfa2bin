
use std::fmt::{Debug};
use std::fs::File;
use std::io::{Write, BufWriter};
use crate::helper::binary2dec_bed;


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

    /// Write bed file
    pub fn write_bed(&self, out_prefix: &str, t: &str){
        //hexdump -C test.bin
        // xxd -b file
        // xxd file
        // SNP: 00000001 , 0
        // IND: 00000000, 1
        let _dummy_vec1: Vec<u8> = vec![108, 27, 0];

        let h = self.copy(1);
        println!("{} {}", h.len(), h[0].len());
        let mut buff: Vec<u8> = vec![108, 27, 1];
        let h2 = transpose(h);
        println!("{} {}", h2.len(), h2[0].len());
        for x in h2.iter(){
            let j: Vec<&[bool]> = x.chunks(4).collect();
            for x in j{
                buff.push(binary2dec_bed(x));
            }
            //println!("Number of bytes {}", buff.len());
            //println!("x {}", x.len());
        }


        let mut file = File::create([out_prefix, t, "bed"].join(".")).expect("Not able to write ");
        file.write_all(&buff).expect("Not able to write ")

    }



}

fn transpose<T>(v: Vec<Vec<T>>) -> Vec<Vec<T>>
    where
        T: Clone,
{
    assert!(!v.is_empty());
    (0..v[0].len())
        .map(|i| v.iter().map(|inner| inner[i].clone()).collect::<Vec<T>>())
        .collect()
}














