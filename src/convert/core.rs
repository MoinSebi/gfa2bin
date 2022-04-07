

use crate::convert::matrix::Matrix;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::io::{Write, BufWriter, BufReader, BufRead};
use std::fs::File;
use gfaR_wrapper::{GraphWrapper, NGfa};
use crate::helper::{binary2dec_bed, binary2dec_bed2, trans2, trans5};
use packing_lib::reader::{ReaderU16, wrapper_bool, ReaderBit, wrapper_u16, get_file_as_byte_vec};
use bimap::{BiMap};
use std::slice::Chunks;
use std::mem;
use bitvec::prelude::*;
use log::info;


#[derive(Debug, Clone, Eq, PartialEq)]
/// Core data structure, which includes ever
pub struct MatrixWrapper {
    pub matrix: Matrix,
    pub column_name: HashMap<u32, String>,

    pub matrix_bin: Vec<bitvec::vec::BitVec>//Vec<Vec<bool>>

}

impl MatrixWrapper {

    /// Dummy initialization
    pub fn new() -> Self {
        let matrx = Matrix::new();
        let col: HashMap<u32, String> = HashMap::new();
        let mut bv: Vec<bitvec::vec::BitVec> = Vec::new();
        Self {
            matrix: matrx,
            column_name: col,
            matrix_bin: bv,
        }
    }


    //--------------------------------------------------------------------------------
    // Modification

    /// Removing genomes
    /// File name includes genome names which should be kept
    pub fn remove_genomes(& mut self, filename: &str){
        let file = File::open(filename).expect("ERROR: CAN NOT READ FILE\n");
        let reader = BufReader::new(file);
        let mut file_genomes: HashSet<String> = HashSet::new();
        let mut bv: Vec<bitvec::vec::BitVec> = Vec::new();

        for line in reader.lines() {
            let l = line.unwrap();
            file_genomes.insert(l);
        }
        let names = file_genomes.clone();


        let names_r:HashSet<String> = self.column_name.values().cloned().collect();
        let kk: HashSet<_> = names_r.difference(&names).collect();


        let mut to_remove = Vec::new();
        let mut rr = Vec::new();
        for x in self.column_name.iter(){
            if names.contains(x.1){
                to_remove.push(x.0.clone());
            }
        }

        for x in self.column_name.iter(){
            if kk.contains(x.1){
                rr.push(x.0.clone());
            }
        }

        info!("{:?}", to_remove);
        info!("this is k {:?}", rr);
        for x in rr.iter(){
            self.column_name.remove(&x);
        }

        rr.sort();
        info!("{:?}", rr);
        for (index, x) in rr.iter().enumerate(){
            self.matrix.matrix_core.remove((x.clone() as usize) - index);
        }
        info!("{:?}", rr);
    }


    /// Split bin matrix into multiple ones
    /// For smaller data and faster read of GEMMA
    pub fn split_bin(&self, number: usize) -> Chunks<bitvec::vec::BitVec>{
        let size = self.matrix_bin.len()/number;
        let mut h = Vec::new();
        let mut tnumb = 0;
        for x in 0..number{
            h.push(&self.matrix_bin[tnumb..tnumb+size]);
        }

        let j = self.matrix_bin.chunks(10);
        j

    }

    /// Split matrix matrix into multiple ones
    /// For smaller data and faster read of GEMMA
    pub fn split_matrix(&self, number: usize) -> Vec<&[Vec<u32>]>{
        let size = self.matrix.matrix_core.len()/number;
        let mut h = Vec::new();
        let mut tnumb = 0;
        for x in 0..number{
            h.push(&self.matrix.matrix_core[tnumb..tnumb+size]);
        }
        let j = self.matrix.matrix_core.chunks(10);
        h

    }

    /// Make a binary matrix (matrix_bin) with a threshold
    /// This is needed for bed output
    pub fn make_binary(& mut self, thresh: u32){
        self.matrix_bin = self.matrix.copy(thresh);
    }


    /// Write BED output
    /// https://zzz.bwh.harvard.edu/plink/binary.shtml
    pub fn write_bed(&self, out_prefix: &str, t: &str){
        //hexdump -C test.bin
        // xxd -b file
        // xxd file
        // SNP: 00000001 , 0
        // IND: 00000000, 1

        let mut buff: Vec<u8> = vec![108, 27, 1];
        // Make SNP Vector
        let h2 = trans5( &self.matrix_bin);
        for x in h2.iter(){
            let j  = x.chunks(4);
            for x in j{
                buff.push(binary2dec_bed2(x));
            }
            //info!("Number of bytes {}", buff.len());
            //info!("x {}", x.len());
        }

        info!("{}", [out_prefix, t, "bed"].join("."));
        let mut file = File::create([out_prefix, t, "bed"].join(".")).expect("Not able to write ");
        file.write_all(&buff).expect("Not able to write ")

    }

    pub fn wrapper_write_bed_split(&self, out_prefix: &str, t: &str){
        let chunks = self.matrix_bin.chunks(10);
        for chunk in chunks{

        }
    }

    /// Filter binary matrix
    /// TODO
    /// Remove transpose and
    /// A lot of transposing -> slow
    pub fn filter(&self) -> Vec<u32>{
        info!("{} {}", self.matrix_bin.len(), self.matrix_bin[0].len());
        let k:Vec<bitvec::vec::BitVec>= trans5(&self.matrix_bin);
        let mut k2 = Vec::new();

        let mut to_remove: Vec<u32> = Vec::new();

        for (i, x) in k.iter().enumerate(){
            let mut sum = 0;
            for y in x.iter() {
                if *y == true {
                    sum += 1;
                }
            }
            if (sum as usize != x.len()) | (sum == 0) {
                k2.push(x.clone());
            } else {
                info!("{} {}", sum, x.len());
                to_remove.push(i as u32);
            }


        }
        info!("Boring SNPs {}", to_remove.len());
        let k3 = trans5(&k2);

        info!("Before {}  After {}", k[0].len(), k3[0].len());
        return to_remove;
    }


    pub fn reduce_combinations_test(& mut self) -> (Vec<usize>, Vec<usize>){
        let mut hm: BiMap<_,_> = BiMap::new();
        let mut h1: Vec<usize> = Vec::new();
        let mut h2: Vec<usize> = Vec::new();
        info!("Starting size2 {}", self.matrix_bin[0].len());

        let mut count = 0;
        for x in 0..self.matrix_bin[0].len(){
            let mut u: bitvec::vec::BitVec = bitvec::vec::BitVec::new();
            for y in 0..self.matrix_bin.len(){
                u.push(self.matrix_bin[y][x]);
            }
            if ! hm.contains_left(&u) {
                hm.insert(u, count);
                h1.push(x);
                h2.push(count);
                count += 1;
            } else {

                h1.push(x);
                h2.push(hm.get_by_left(&u).unwrap().clone());
            }

        }
        let mut h : Vec<bitvec::vec::BitVec> = Vec::new();
        for x in 0..hm.iter().len(){
            h.push(hm.get_by_right(&x).unwrap().clone());
        }
        info!("Starting size2 {}", self.matrix_bin[0].len());
        self.matrix_bin = trans5(&h);


        info!("djsakdhsja {}", self.matrix_bin[0].len());

        (h1, h2)
    }

    #[allow(dead_code)]
    /// Writing bim file
    /// Information: https://www.cog-genomics.org/plink/1.9/formats#bim
    pub fn write_bim(&self, out_prefix: &str, t: &str) {


        let f = File::create([out_prefix, t,  "bim"].join(".")).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        for x in 0..self.matrix_bin[0].len(){
            write!(f, "{}\t{}\t{}\t{}\t{}\t{}\n", "graph", ".", 0, x, "A", "T").expect("Not able to write ");
        }

    }


    /// Reduce binary shit
    pub fn reduce_combinations(& mut self) -> (Vec<usize>, Vec<usize>){
        // Meta
        // h1,h2 -> meta
        // hm -> BiMap (vec -> usize)
        let mut hm: BiMap<_,_> = BiMap::new();
        let mut h1: Vec<usize> = Vec::new();
        let mut h2: Vec<usize> = Vec::new();


        // Make SNPs Vector
        info!("reduce it djsakdhsja {}", self.matrix_bin.len());
        let k: Vec<bitvec::vec::BitVec>= trans5(&self.matrix_bin);
        info!("reduce it djsakdhsja2 {}", self.matrix_bin.len());
        info!("Starting size {}", k.len());

        let mut count = 0;
        // Iterate over SNPs
        for (index, x) in k.iter().enumerate() {
            if ! hm.contains_left(x) {
                hm.insert(x.clone(), count);
                h1.push(index);
                h2.push(count);
                count += 1;
            } else {

                h1.push(index);
                h2.push(hm.get_by_left(x).unwrap().clone());
            }
        }

        // Make a Vector (from 1 - len(x))
        let mut h : Vec<bitvec::vec::BitVec> = Vec::new();
        for x in 0..hm.iter().len(){
            h.push(hm.get_by_right(&x).unwrap().clone());
        }

        // Make Acc vector
        info!("Reduced size {:?}", h.len());
        let h_convert = trans5(&h);
        self.matrix_bin = h_convert;
        (h1,h2)
    }


    /// Write the names - helper function
    pub fn write_names(&self, out_prefix: &str) {
        let f = File::create([out_prefix,  "bim_names"].join(".")).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        for x in self.column_name.iter(){
            write!(f, "{}\n", x.1).expect("Can not write file");
        }
    }



}



//---------------------------------------------------------------------------------------------------
// This is for writing (may move later)
#[allow(dead_code)]
/// Writing bim file
/// Information: https://www.cog-genomics.org/plink/1.9/formats#bim
pub fn write_bim<T>(ll: &BiMap<T, usize>, out_prefix: &str, t: &str)
    where
        T: Debug + std::hash::Hash + std::cmp::Eq + Ord
{


    let f = File::create([out_prefix, t,  "bim"].join(".")).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    for x in 0..ll.len(){
        write!(f, "{}\t{}\t{}\t{}\t{}\t{}\n", "graph", ".", 0, x, "A", "T").expect("Not able to write ");
    }

}

/// Writing bim helper
/// Index -> Feature
/// Index because wrongly removed
pub fn write_bimhelper<T>(ll: &BiMap<T, usize>, out_prefix: &str, t: &str)
    where
        T: Debug + std::hash::Hash + std::cmp::Eq + Ord
{


    let f = File::create([out_prefix, t,  "bimhelper"].join(".")).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    let o: Vec<&usize> = ll.right_values().collect();
    for x in 0..ll.right_values().len(){
        write!(f, "{}\t{:?}\n", x, ll.get_by_right(&x).unwrap()).expect("Not able to write ");
    }

}



//------------------------------------------------------------------------------------------------------------------------
// Modification

/// Remove entries from bimap
pub fn remove_bimap<T>(bm: & mut BiMap<T, usize>, v: Vec<u32>)
    where
        T:  Debug + std::hash::Hash + std::cmp::Eq
{

    for x in v.iter(){
        bm.remove_by_right(&(*x as usize));
    }

}

