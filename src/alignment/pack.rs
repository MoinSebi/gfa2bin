use bitvec::order::Msb0;
use bitvec::prelude::BitVec;
use log::info;
use packing_lib::reader::{get_file_as_byte_vec, ReaderBit, ReaderU16, wrapper_bool, wrapper_u16};
use crate::core::core::MatrixWrapper;
use crate::graph::parser::Haplotype;

#[allow(dead_code)]
/// Matrix constructor from mappings (bits)
pub fn matrix_pack_bit(filename: &str, matrix_w: & mut MatrixWrapper, h2: &mut Vec<u32>) {
    let buf: Vec<u8> = get_file_as_byte_vec(filename);
    let data: Vec<ReaderBit> = wrapper_bool(&buf);
    for (index, reader_bit) in data.iter().enumerate(){
        //println!("{}", x.name);
        matrix_w.core_names.insert(index as u32, reader_bit.name.clone());
        //println!("{}", k.len());
        let k = reader_bit.data.clone();
        matrix_w.matrix_bin.push(k);
    }
    let length = matrix_w.matrix_bin[0].len();
    h2.extend(0..length as u32);
}

pub fn matrix_pack_bit_v2(filename: &str, matrix_w: & mut MatrixWrapper, h2: &mut Vec<u32>) {
    let buf: Vec<u8> = get_file_as_byte_vec(filename);
    let k: Vec<ReaderBit> = wrapper_bool(&buf);
    let hap = Haplotype::from_string_vec(k.iter().map(|x| x.name.clone()).collect(), "#");
    let ll = k[0].data.len().clone();
    matrix_w.matrix_bin = vec![BitVec::<u8, Msb0>::repeat(false, hap.haplotype.len() * 2); ll];
    for (i2, x) in hap.haplotype2.iter().enumerate(){
        matrix_w.bin_names.push(x.0.clone());
        for (i, y) in k[x.1[0]].data.iter().enumerate(){
            if y == &true {
                matrix_w.matrix_bin[i].get_mut(i2*2).unwrap().set(true);
            }
        }
        for (i, y) in k[x.1[1]].data.iter().enumerate(){
            if y == &true {
                matrix_w.matrix_bin[i].get_mut(i2*2+1).unwrap().set(true);
            }
        }
    }


    h2.extend(1..(matrix_w.matrix_core.len() +1) as u32);
}


/// Matrix constructor from mappings (u16)
pub fn matrix_pack_u16_v2(filename: &str, matrix_w: & mut MatrixWrapper, h2: & mut Vec<u32>) {
    let g: Vec<u8> = get_file_as_byte_vec(filename);
    let k: Vec<ReaderU16> = wrapper_u16(&g);
    let ll = k[0].data.len().clone();
    matrix_w.matrix_core = vec![vec![0; k.len()]; ll];
    for (i, x) in k.iter().enumerate(){
        matrix_w.names.push(x.name.clone());
        for (i2, y) in x.data.iter().enumerate(){
            matrix_w.matrix_core[i2][i] = *y as u32;
        }
    }

    h2.extend(1..(matrix_w.matrix_core.len() +1) as u32);
}

#[allow(dead_code)]
pub fn matrix_pack_u16(filename: &str, matrix_w: & mut MatrixWrapper, h2: & mut Vec<u32>) {
    let g: Vec<u8> = get_file_as_byte_vec(filename);
    let k: Vec<ReaderU16> = wrapper_u16(&g);

    for (index, readeru16) in k.iter().enumerate(){
        //println!("{}", x.name);
        matrix_w.core_names.insert(index as u32, readeru16.name.clone());
        // First map function use!
        let u: Vec<u32> = readeru16.data.clone().iter().map(|f| f.clone() as u32).collect();
        matrix_w.matrix_core.push(u);
    }
    info!("Make BIMAP");
    let length = matrix_w.matrix_core[0].len();

    h2.extend(0..length as u32);
}



