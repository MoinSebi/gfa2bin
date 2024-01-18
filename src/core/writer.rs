use std::fmt::{Debug};
use std::fs::File;
use std::io::{Write, BufWriter};
use bitvec::order::Msb0;
use bitvec::vec::BitVec;
use crate::helper::{bitvec_to_u8};













pub fn write_bed2(sel: &[BitVec<u8, Msb0>], out_prefix: &str, feature: &str, number: usize, total_len: usize){
    //hexdump -C test.bin
    // xxd -b file
    // xxd file
    // SNP: 00000001 , 0
    // IND: 00000000, 1

    let mut buff: Vec<u8> = vec![108, 27, 1];
    // Make SNP Vector
    for x in sel.iter(){
        let j  = x.chunks(8);
        for x in j{
            buff.push(bitvec_to_u8(x));
        }
        //info!("Number of bytes {}", buff.len());
        //info!("x {}", x.len());
    }

    let mut output = [out_prefix,  feature, &number.to_string(), "bed"].join(".");
    if total_len == 1{
        output = [out_prefix,  feature, "bed"].join(".");
    }
    let mut file = File::create(output).expect("Not able to write ");
    file.write_all(&buff).expect("Not able to write ")

}



/// Write a bim file to a file
///
/// Based on the number of path/samples as input
/// Based on nodes
pub fn write_bim_nodes(names: &[usize], out_prefix: &str, feature: &str, number: usize, totallen: usize){
    let mut output = [out_prefix,  feature, &number.to_string(), "bim"].join(".");
    if totallen == 1{
        output = [out_prefix,  feature, "bim"].join(".");
    }
    let f = File::create(output).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    for x in names.iter() {
        write!(f, "{}\t{}\t{}\t{}\t{}\t{}\n", "graph", ".", 0, x, "A", "T").expect("Can not write file");
    }

}


/// Write a bim file to a file
///
/// Based on the number of path/samples as input
/// Based on node_dir
pub fn write_bim_dirnode(names: &[(usize, bool)], out_prefix: &str, feature: &str, number: &usize, iter: usize){
    let mut output = [out_prefix,  feature, &number.to_string(), "bim"].join(".");
    if iter == 1{
        output = [out_prefix,  feature, "bim"].join(".");
    }
    let f = File::create(output).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    for x in names.iter() {
        write!(f, "{}\t{}\t{}\t{}\t{}\t{}\n", "graph", ".", 0, x.0.to_string() + &(if x.1 { "+" } else { "-" }).to_string() , "A", "T").expect("Can not write file");
    }

}

/// Write a bim file to a file
///
/// Based on the number of path/samples as input
/// Based on edges
pub fn write_bim_edges(names: &[(u32, bool, u32, bool)], out_prefix: &str, feature: &str, number: &usize, iter: usize){
    let mut output = [out_prefix,  feature, &number.to_string(), "bim"].join(".");
    if iter == 1{
        output = [out_prefix,  feature, "bim"].join(".");
    }
    let f = File::create(output).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    for x in names.iter() {
        write!(f, "{}\t{}\t{}\t{}\t{}\t{}\n", "graph", ".", 0, x.0.to_string() + &(if x.1 { "+" } else { "-" }).to_string() + &x.2.to_string() + &(if x.3 { "+" } else { "-" }).to_string(), "A", "T").expect("Can not write file");
    }

}

