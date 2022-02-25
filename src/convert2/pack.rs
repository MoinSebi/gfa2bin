use bimap::BiMap;
use packing_lib::reader::{get_file_as_byte_vec, ReaderBit, ReaderU16, wrapper_bool, wrapper_u16};
use crate::MatrixWrapper;

/// Matrix constructor from mappings (u16)
pub fn matrix_pack_bit(filename: &str, mw: & mut MatrixWrapper, h2: & mut BiMap<u32, usize>) {
    let g: Vec<u8> = get_file_as_byte_vec(filename);
    let k: Vec<ReaderBit> = wrapper_bool(&g);

    for (i,x) in k.iter().enumerate(){
        //println!("{}", x.name);
        mw.column_name.insert(i as u32,x.name.clone());
        //println!("{}", k.len());
        mw.matrix_bin.push(x.cc.clone());
    }
    eprintln!("Make BIMAP");
    eprintln!("size {}", mw.matrix_bin[0].len());
    for x in 0..mw.matrix_bin[0].len(){
        h2.insert(x as u32, x);
    }
}

pub fn makeBIMAP(maxlen: usize){
    let mut j1: hashbrown::HashMap<u32, usize> = hashbrown::HashMap::new();
    let mut j2: hashbrown::HashMap<usize, u32> = hashbrown::HashMap::new();
    for x in 0..maxlen{
        j1.insert(x as u32, x.clone());
        j2.insert(x, x as u32);
    }

}


/// Matrix constructor from mappings (u16)
pub fn matrix_pack_u16(filename: &str, mw: & mut MatrixWrapper, h2: & mut BiMap<u32, usize>) {
    let g: Vec<u8> = get_file_as_byte_vec(filename);
    let k: Vec<ReaderU16> = wrapper_u16(&g);

    for (i,x) in k.iter().enumerate(){
        //println!("{}", x.name);
        mw.column_name.insert(i as u32,x.name.clone());
        // First map function use!
        let u: Vec<u32> = x.cc.clone().iter().map(|f| f.clone() as u32).collect();
        mw.matrix.matrix_core.push(u);
    }
    eprintln!("Make BIMAP");
    for x in 0..mw.matrix.matrix_core[0].len(){
        h2.insert(x as u32, x);
    }
}
