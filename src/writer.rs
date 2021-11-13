use std::fs::File;
use std::io::{Write, BufWriter};
use bimap::BiMap;
use std::fmt::Debug;

/// Write the names - helper function
pub fn write_reduce(h1: &Vec<usize>, h2:  &Vec<usize>, out_prefix: &str) {
    let f = File::create([out_prefix,  "reduce"].join(".")).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    for x in 0..h1.len(){
        write!(f, "{}\t{}\n", h1[x], h2[x]).expect("Can not write file");
    }
}

pub fn write_bimap<T>(bm: &BiMap<T, usize>)
    where
        T:  Debug + std::hash::Hash + std::cmp::Eq
{
    let f = File::create("bimbim").expect("Unable to create file");
    let mut f = BufWriter::new(f);
    for (k1,k2) in bm.iter(){
        write!(f, "{:?}\t{}\n", k1, k2).expect("Can not write file");
    }
}