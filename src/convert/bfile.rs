use std::fs;
use std::fs::File;
use std::io::{BufReader, BufRead, Read};
use bitvec::order::Msb0;
use bitvec::prelude::BitVec;
use crate::MatrixWrapper;


/// Wrapper for reading all entries of a bfile
/// matrix_bin -> bim, bam
/// column_name
/// helper
pub fn bfile_wrapper(filename: &str, matrix: &mut MatrixWrapper, names: &mut Vec<String>) {
    matrix.transposed = true;
    read_fam(&format!("{}{}", filename, ".fam"), matrix);
    read_bim(&format!("{}{}", filename, ".bim"), names);

    read_bed(&format!("{}{}", filename, ".bed"), matrix,  names.len());
}

/// Read fam file and extract only the names
/// Add it to the matrix
pub fn read_fam(filename: &str, matrix: &mut MatrixWrapper){
    let data = fs::read_to_string(filename).expect("Unable to read file");
    let dd: Vec<_> = data.split("\n").collect();
    println!("{}", dd.len());
    let dd2: Vec<_> = dd.into_iter().map(| e| e.to_string()).collect();
    for (i, x) in dd2.into_iter().enumerate(){
        let gg: Vec<_> = x.split("\t").collect();
        matrix.core_names.insert(i as u32, gg[0].to_string());
    }
}

// test, helper (chr_pos) (string - might be big lol)
pub fn read_bim(filename: &str, snp_names: &mut Vec<String>){
    let file = BufReader::new(File::open(filename).expect("Unable to open file"));


    for x in file.lines() {
        let f: Vec<String> = x.unwrap().split("\t").map(| s | s.to_string()).collect();
        let ff = f[0].clone();
        let ff2 = f[3].clone();

        snp_names.push(format!("{}_{}", ff.to_string(), ff2));
    }
    println!("dasd{}", snp_names.len());
}

pub fn read_bed(filename: &str, matrix_w: & mut MatrixWrapper, numbsnp: usize) {

    let mut file = File::open(filename).expect("no file found");
    let metadata = fs::metadata(&filename).expect("unable to read metadata");
    let mut buffer = vec![0; metadata.len() as usize];

    file.read_exact(&mut buffer).expect("buffer overflow");
    // cut off first 3 bytes
    buffer = buffer[3..].to_vec();

    // do i need this
    let mut num = matrix_w.core_names.len()/4;
    if (matrix_w.core_names.len() % 4) > 0{
        num += 1;
    }
    // Each chunk
    let chunks = buffer.chunks(num);

    println!("{}", numbsnp);
    for chunk in chunks.into_iter() {
        let bv: BitVec<u8, Msb0> = BitVec::from_slice(&chunk[..]);
        let mut dd: BitVec<u8, Msb0> = BitVec::new();

        for (i, x) in bv.iter().step_by(2).enumerate(){
            if i < num as usize {
                dd.push(*x);
            }

        }
        matrix_w.matrix_bin.push(dd);
    }
}

