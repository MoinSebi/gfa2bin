use std::fmt::Debug;
use std::fs::File;
use std::io::{Write, BufWriter, BufReader, BufRead};

use bimap::BiMap;
use crate::MatrixWrapper;


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
    for x in 0..ll.right_values().len(){
        write!(f, "{}\t{:?}\n", x, ll.get_by_right(&x).unwrap()).expect("Not able to write ");
    }

}





/// Writing file wrapper
/// BED + BIM + GENOME ORDER
pub fn write_matrix(se: & mut MatrixWrapper, what: &str, out_prefix: &str, t: &str){
    if (what == "bed") | (what == "all"){
        if se.matrix_bin.is_empty(){
            se.matrix_bin = se.matrix.copy(1);
        };
        se.write_bed(out_prefix, t);
        se.write_bim(out_prefix, "gfa2bin");

        write_genome_order(se, out_prefix);
    }
    if (what == "bimbam") | (what == "all"){
        //se.matrix.write_bimbam(out_prefix, t);
    }
}



/// Write the order of genomes
pub fn write_genome_order(se: & mut MatrixWrapper, out_prefix: &str){
    let f = File::create([out_prefix,  "bim_names"].join(".")).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    for x in 0..se.column_name.len(){
        write!(f, "{}\n", se.column_name.get(&(x as u32)).unwrap()).expect("Can not write file");
    }
}


use crate::helper::{trans2, binary2dec_bed, trans3};
use std::slice::Chunks;

/// Write the names - helper function
pub fn write_reduce(h1: &Vec<usize>, h2:  &Vec<usize>, out_prefix: &str, t: &str) {
    let f = File::create([out_prefix,  t, "reduce"].join(".")).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    for x in 0..h1.len(){
        write!(f, "{}\t{}\n", h1[x], h2[x]).expect("Can not write file");
    }
}

/// Write BIMAP (index)
/// This function may be redundant in the future
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

/// Write BIMAP (but you can split)
pub fn write_bimap2<T>(bm: &BiMap<T, usize>, till: usize, number: usize)
    where
        T:  Debug + std::hash::Hash + std::cmp::Eq
{
    let f = File::create("bimbim").expect("Unable to create file");
    let mut f = BufWriter::new(f);
    for (index, (k1,k2)) in bm.iter().enumerate(){
        if index > till{
            break
        }
        write!(f, "{:?}\t{}\n", k1, k2).expect("Can not write file");
    }
}


/// For multiple bed files
/// Splitting
pub fn write_bed_split(data: &[bitvec::vec::BitVec], out_prefix: &str, t: &str){
    //hexdump -C test.bin
    // xxd -b file
    // xxd file
    // SNP: 00000001 , 0
    // IND: 00000000, 1

    let mut buff: Vec<u8> = vec![108, 27, 1];
    // Make SNP Vector
    let h2 = trans3( &data);
    for x in h2.iter(){
        let j: Vec<&[bool]> = x.chunks(4).collect();
        for x in j{
            buff.push(binary2dec_bed(x));
        }
        //println!("Number of bytes {}", buff.len());
        //println!("x {}", x.len());
    }

    println!("{}", [out_prefix, t, "bed"].join("."));
    let mut file = File::create([out_prefix, t, "bed"].join(".")).expect("Not able to write ");
    file.write_all(&buff).expect("Not able to write ")

}

