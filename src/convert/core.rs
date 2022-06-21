

use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::io::{Write, BufWriter, BufReader, BufRead};
use std::fs::File;
use bimap::{BiMap};
use std::slice::Chunks;
use bitvec::prelude::*;
use log::info;
use crate::helper::{binary2dec_bed2, trans5};


#[derive(Debug, Clone, Eq, PartialEq)]
/// Core data structure, which includes ever
pub struct MatrixWrapper {
    pub shape: (usize, usize),
    pub matrix_core: Vec<Vec<u32>>,
    pub column_name: HashMap<u32, String>,
    pub matrix_bin: Vec<BitVec<u8, Msb0>>//Vec<Vec<bool>>

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
            column_name: col,
            matrix_bin: bv,
        }
    }


    //--------------------------------------------------------------------------------
    pub fn make_binary(& mut self, number: u32){
        let mut new_matrix: Vec<BitVec<u8, Msb0>> = Vec::new();
        for x in self.matrix_core.iter(){
            let mut new_vec: BitVec<u8, Msb0> = BitVec::new();
            for y in x.iter(){
                new_vec.push(y.clone() >= number);

            }
            new_matrix.push(new_vec);
        }
        self.matrix_bin = new_matrix;
    }


    //--------------------------------------------------------------------------------
    // Modification

    /// Removing genomes
    /// File name includes genome names which should be kept
    pub fn remove_genomes(& mut self, filename: &str){
        let file = File::open(filename).expect("ERROR: CAN NOT READ FILE\n");
        let reader = BufReader::new(file);
        let mut file_genomes: HashSet<String> = HashSet::new();

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
            self.matrix_core.remove((x.clone() as usize) - index);
        }
        info!("{:?}", rr);
    }


    /// Split bin matrix into multiple ones
    /// For smaller data and faster read of GEMMA
    pub fn split_bin(&self, number: usize) -> Chunks<BitVec<u8, Msb0>>{
        let size = ((self.matrix_bin.len() as f64)/(number as f64)).ceil() as usize;
        let j = self.matrix_bin.chunks(size);
        j

    }

    #[allow(dead_code)]
    #[allow(dead_code)]
    /// Split matrix matrix into multiple ones
    /// For smaller data and faster read of GEMMA
    pub fn split_matrix(&self, number: usize) -> Vec<&[Vec<u32>]>{
        let size = self.matrix_core.len()/number;
        let mut h = Vec::new();
        let mut tnumb = 0;
        for _x in 0..number{
            h.push(&self.matrix_core[tnumb..tnumb+size]);
            tnumb += tnumb+size;
        }
        h

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
        let h2 = &self.matrix_bin;
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


    /// Filter the (all-shared) entries
    ///
    /// Returns: Vector of index of entries which are shared by all.
    pub fn filter_shared(&self) -> Vec<u32>{
        info!("{} {}", self.matrix_bin.len(), self.matrix_bin[0].len());
        let mut k2 = Vec::new();

        let mut to_remove: Vec<u32> = Vec::new();

        for (i, x) in self.matrix_bin.iter().enumerate(){
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

        return to_remove;
    }


    pub fn reduce_combinations_test(& mut self) -> (Vec<usize>, Vec<usize>){
        let mut hm: BiMap<_,_> = BiMap::new();
        let mut h1: Vec<usize> = Vec::new();
        let mut h2: Vec<usize> = Vec::new();
        info!("Starting size2 {}", self.matrix_bin[0].len());

        let mut count = 0;
        for x in 0..self.matrix_bin[0].len(){
            let mut u: BitVec<u8, Msb0> = BitVec::new();
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
        let mut h : Vec<BitVec<u8, Msb0>> = Vec::new();
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
        for x in 0..self.matrix_bin.len(){
            write!(f, "{}\t{}\t{}\t{}\t{}\t{}\n", "graph", ".", 0, x, "A", "T").expect("Not able to write ");
        }

    }

    #[allow(dead_code)]
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
        let k: Vec<BitVec<u8, Msb0>>= trans5(&self.matrix_bin);
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
        let mut h : Vec<BitVec<u8, Msb0>> = Vec::new();
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

#[allow(dead_code)]
/// Writing bim helper
/// Index -> Feature
/// Index because wrongly removed
pub fn write_bimhelper<T>(ll: &BiMap<T, usize>, out_prefix: &str, t: &str)
    where
        T: Debug + std::hash::Hash + std::cmp::Eq + Ord
{


    let f = File::create([out_prefix, t,  "bimhelper"].join(".")).expect("Unable to create file");
    let mut f = BufWriter::new(f);
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

    println!("{}", bm.len());
    for x in v.iter(){
        bm.remove_by_right(&(*x as usize));
    }
    println!("{}", bm.len());


}

pub fn write_bim2(start: usize, end: usize, out_prefix: &str, t: &str) {


    let f = File::create([out_prefix, t,  "bim"].join(".")).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    for x in start..end{
        write!(f, "{}\t{}\t{}\t{}\t{}\t{}\n", "graph", ".", 0, x, "A", "T").expect("Not able to write ");
    }

}

