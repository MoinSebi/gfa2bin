use bitvec::order::{Msb0};
use bitvec::prelude::BitVec;
use bitvec::slice::BitSlice;
use byteorder::{BigEndian, ByteOrder};
use log::info;
use packing_lib::reader::get_file_as_byte_vec;


#[allow(dead_code)]
pub fn binary2dec_bed(vecc: &[bool]) -> u8{
    let mut result: u8 = 0;
    let mut count = 0;
    for x in vecc.iter(){
        let t: u8 = 2;
        result += (t.pow(count as u32)) * (*x as u8);
        count += 1;
        result += (t.pow(count as u32)) * (*x as u8);
        count += 1;
    }
    result
}




pub fn bitvec_to_u8(bitvec: &BitSlice<u8, Msb0>) -> u8 {
    let mut result: u8 = 0;

    // Calculate the number of zero bits to pad on the left
    let padding_bits = 8 - bitvec.len();

    for (i, bit) in bitvec.iter().enumerate() {
        if *bit {
            result |= 1 << (padding_bits + i);
        }
    }

    result
}

pub fn make_dir_name(maxval: &usize) -> Vec<(usize, bool)>{
    let mut f = Vec::new();
    for x in 0..*maxval{
        f.push((x, false));
        f.push((x, true));
    }
    f
}

pub fn is_all_zeros(bitvector: &BitVec<u8, Msb0>) -> bool {
    for byte in bitvector.iter() {
        if byte.as_ref() != &false {
            return false; // If any byte is not zero, the bitvector is not all zeros
        }
    }
    true // If no non-zero byte is found, the bitvector is all zeros
}

pub fn is_all_ones(bitvector: &BitVec<u8, Msb0>) -> bool {
    for byte in bitvector.iter() {
        if byte.as_ref() != &true {
            return false; // If any byte is not zero, the bitvector is not all zeros
        }
    }
    true // If no non-zero byte is found, the bitvector is all zeros
}

pub fn custom_retain_two_vectors<T, F, S>(vec1: &mut Vec<T>, vec2: &mut Vec<S>, condition: F)
    where
        F: Fn(&T) -> bool,
{
    assert_eq!(vec1.len(), vec2.len(), "Vectors must have the same length");

    let mut write_index = 0;

    for read_index in 0..vec1.len() {
        if condition(&vec1[read_index]) {
            // Retain elements that satisfy the condition
            if write_index != read_index {
                // Move the elements to their new positions
                vec1.swap(write_index, read_index);
                vec2.swap(write_index, read_index);
            }

            write_index += 1;
        }
    }

    // Truncate both vectors to their new length
    vec1.truncate(write_index);
    vec2.truncate(write_index);
}

#[allow(dead_code)]
/// Transposing a 2D vector
///
/// Creating a totally new 2D vector
pub fn transpose_generic<T>(v: &Vec<Vec<T>>) -> Vec<Vec<T>>
    where
        T: Clone,
{
    info!("Transposing generic");
    assert!(!v.is_empty());
    (0..v[0].len())
        .map(|i| v.iter().map(|inner| inner[i].clone()).collect::<Vec<T>>())
        .collect()
}

#[allow(dead_code)]
/// Transposing a 2D Vec<BitVec>
///
/// Creating a totally new vector
pub fn transpose_bitvec(v: &Vec<BitVec<u8, Msb0>>) -> Vec<BitVec<u8, Msb0>>{
    info!("Transposing bitvec");
    let mut o: Vec<BitVec<u8, Msb0>> = Vec::new();
    for x in 0..v[0].len(){
        let mut o2: BitVec<u8, Msb0> = BitVec::new();
        for y in 0..v.len(){
            o2.push(v[y][x].clone());
        }
        o.push(o2);
    }
    return o;
}

pub fn get_thresh(filename: &str) -> u16{
    let size = BigEndian::read_u16(&mut & get_file_as_byte_vec(filename)[7..9]);
    size
}
