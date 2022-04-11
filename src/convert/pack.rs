use bimap::BiMap;
use log::{info};
use packing_lib::reader::{get_file_as_byte_vec, ReaderBit, ReaderU16, wrapper_bool, wrapper_u16};
use crate::MatrixWrapper;

/// Matrix constructor from mappings (bits)
pub fn matrix_pack_bit(filename: &str, matrix_w: & mut MatrixWrapper, h2: & mut BiMap<u32, usize>) {
    let buf: Vec<u8> = get_file_as_byte_vec(filename);
    let data: Vec<ReaderBit> = wrapper_bool(&buf);
    for (index, reader_bit) in data.iter().enumerate(){
        //println!("{}", x.name);
        matrix_w.column_name.insert(index as u32, reader_bit.name.clone());
        //println!("{}", k.len());
        let k = reader_bit.data.clone();
        matrix_w.matrix_bin.push(k);
    }
    for x in 0..matrix_w.matrix_bin[0].len(){
        h2.insert(x as u32, x);
    }
}

#[allow(dead_code)]
pub fn make_bimap(maxlen: usize){
    let mut j1: hashbrown::HashMap<u32, usize> = hashbrown::HashMap::new();
    let mut j2: hashbrown::HashMap<usize, u32> = hashbrown::HashMap::new();
    for x in 0..maxlen{
        j1.insert(x as u32, x.clone());
        j2.insert(x, x as u32);
    }

}


/// Matrix constructor from mappings (u16)
pub fn matrix_pack_u16(filename: &str, matrix_w: & mut MatrixWrapper, h2: & mut BiMap<u32, usize>) {
    let g: Vec<u8> = get_file_as_byte_vec(filename);
    let k: Vec<ReaderU16> = wrapper_u16(&g);

    for (index, readeru16) in k.iter().enumerate(){
        //println!("{}", x.name);
        matrix_w.column_name.insert(index as u32, readeru16.name.clone());
        // First map function use!
        let u: Vec<u32> = readeru16.data.clone().iter().map(|f| f.clone() as u32).collect();
        matrix_w.matrix_core.push(u);
    }
    info!("Make BIMAP");
    for x in 0..matrix_w.matrix_core[0].len(){
        h2.insert(x as u32, x);
    }
}
