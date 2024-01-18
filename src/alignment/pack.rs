use crate::core::core::MatrixWrapper;
use bitvec::order::{Lsb0, Msb0};
use bitvec::prelude::BitVec;
use packing_lib::reader::{get_file_as_byte_vec, wrapper_bool, wrapper_u16, ReaderBit, ReaderU16};

pub fn matrix_pack_bit_v2(filename: &str, matrix_w: &mut MatrixWrapper, h2: &mut Vec<u32>) {
    let buf: Vec<u8> = get_file_as_byte_vec(filename);
    let k: Vec<ReaderBit> = wrapper_bool(&buf);
    let ll = k[0].data.len().clone();

    matrix_w.matrix_bin = vec![BitVec::<u8, Lsb0>::repeat(false, k.len() * 2); ll];
    for (i2, x) in k.iter().enumerate() {
        matrix_w.sample_names.push(x.name.clone());
        for (i, y) in x.data.iter().enumerate() {
            if y == &true {
                matrix_w.matrix_bin[i].get_mut(i2 * 2).unwrap().set(true);
                matrix_w.matrix_bin[i]
                    .get_mut(i2 * 2 + 1)
                    .unwrap()
                    .set(true);
            }
        }
    }

    h2.extend(1..(matrix_w.matrix_core.len() + 1) as u32);
}

/// Matrix constructor from mappings (u16)
pub fn matrix_pack_u16_v2(filename: &str, matrix_w: &mut MatrixWrapper, h2: &mut Vec<u32>) {
    let g: Vec<u8> = get_file_as_byte_vec(filename);
    let k: Vec<ReaderU16> = wrapper_u16(&g);
    let ll = k[0].data.len().clone();
    matrix_w.matrix_core = vec![vec![0; k.len()]; ll];
    for (i, x) in k.iter().enumerate() {
        matrix_w.sample_names.push(x.name.clone());
        for (i2, y) in x.data.iter().enumerate() {
            matrix_w.matrix_core[i2][i] = *y as u32;
        }
    }
    h2.extend(1..(matrix_w.matrix_core.len() + 1) as u32);
}
