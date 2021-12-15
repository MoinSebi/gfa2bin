use std::fs::File;
use std::io::{Write, BufWriter};
use bimap::BiMap;
use std::fmt::Debug;
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
pub fn write_bed_split(data: &[Vec<bool>], out_prefix: &str, t: &str){
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

/// For multiple bed files
/// Splitting
pub fn write_bed_split2(data: &[Vec<bool>], out_prefix: &str, t: &str){
    //hexdump -C test.bin
    // xxd -b file
    // xxd file
    // SNP: 00000001 , 0
    // IND: 00000000, 1

    let mut buff: Vec<u8> = vec![108, 27, 1];
    // Make SNP Vector
    let h2 = trans3( &data);
    let mut o: Vec<Vec<bool>> = Vec::new();
    for x in 0..data[0].len(){
        let mut o2: Vec<bool> = Vec::new();
        for y in 0..v.len(){
            o2.push(v[y][x].clone());
        }
        let j: Vec<&[bool]> = o2.chunks(4).collect();
        for x2 in j{
            buff.push(binary2dec_bed(x2));
        }

    }


    println!("{}", [out_prefix, t, "bed"].join("."));
    let mut file = File::create([out_prefix, t, "bed"].join(".")).expect("Not able to write ");
    file.write_all(&buff).expect("Not able to write ")


}