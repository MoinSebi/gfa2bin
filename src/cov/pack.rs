use crate::core::core::MatrixWrapper;
use crate::core::helper::{index2node_seq, merge_u32_to_u64};
use bitvec::order::Lsb0;
use bitvec::prelude::BitVec;


use packing_lib::core::core::{DataType, PackCompact};
use packing_lib::normalize::convert_helper::Method;

/// # Binary Pack to Binary Matrix
pub fn bin2bin(matrix_w: &mut MatrixWrapper, input: &PackCompact, index: usize) {
    for (i, y) in input.bin_coverage.iter().enumerate() {
        if y == &true {
            matrix_w.matrix_bit[i].get_mut(index * 2).unwrap().set(true);
            matrix_w.matrix_bit[i]
                .get_mut(index * 2 + 1)
                .unwrap()
                .set(true);
        }
    }
}

/// # f32 Pack to Binary Matrix
pub fn f32_to_bin(matrix_w: &mut MatrixWrapper, input: &PackCompact, thresh: f32, index: usize) {
    for (i, y) in input.normalized_coverage.iter().enumerate() {
        if y > &thresh {
            matrix_w.matrix_bit[i].get_mut(index * 2).unwrap().set(true);
            matrix_w.matrix_bit[i]
                .get_mut(index * 2 + 1)
                .unwrap()
                .set(true);
        }
    }
}

/// # f32 Pack to f32 Matrix
pub fn f32_to_f32(matrix_w: &mut MatrixWrapper, input: &PackCompact, thresh: f32, index: usize) {
    for (i, y) in input.normalized_coverage.iter().enumerate() {
        if y > &thresh {
            matrix_w.matrix_f32[i][index] = *y;
        }
    }
}

/// # Read pack wrapper
///
/// - If the pack is plain, parse it
/// - If the pack is compressed, read it
pub fn read_pack_wrapper(is_plain: bool, filename: &str) -> PackCompact {
    if is_plain {
        PackCompact::parse_pack(filename)
    } else {
        PackCompact::read_wrapper(filename)
    }
}

pub fn wrapper_reader123(buffer: &[u8]) -> PackCompact {
    // total length 85 + len
    let (_kind, _include_all, _bin, _method, _relative, _std, _thresh, _bytes, _length, _name) =
        PackCompact::get_meta(buffer);
    if _bin == DataType::TypeU16 {
        PackCompact::read_u16(buffer)
    } else if _bin == DataType::TypeBit {
        PackCompact::read_bin_coverage(buffer)
    } else {
        PackCompact::read_f32(buffer)
    }
}

/// # Init geno names by index
///
/// Use either
/// - `index` if `pc.node_index` is empty
/// - `pc.node_index` if `pc.node_index` is not empty
pub fn init_geno_names(
    mw: &mut MatrixWrapper,
    pc: &mut PackCompact,
    want_node: bool,
    index: &Vec<u32>,
) {
    /// Check the index
    if pc.node_index.is_empty() {
        pc.node_index = index.clone();
    }

    if want_node {
        mw.geno_names = remove_duplicates(&pc.node_index);
    } else {
        mw.geno_names = index2node_seq(&pc.node_index)
    }
}

/// # Init matrix
///
/// Based on the output
pub fn init_matrix(
    mw: &mut MatrixWrapper,
    pc: &mut PackCompact,
    want_node: bool,
    bimbam: bool,
    len1: usize,
) {
    // Create normalized coverage from coverage
    if pc.bin_coverage.is_empty() {
        if pc.normalized_coverage.is_empty() {
            if want_node {
                pc.calc_node_cov();
            } else {
                pc.normalized_coverage = pc.coverage.iter().map(|x| *x as f32).collect();
            }
        } else if want_node {
            pc.calc_node_cov();
        }

        // if bin_coverage is there
    } else {
        mw.matrix_bit = vec![BitVec::<u8, Lsb0>::repeat(false, len1 * 2); pc.bin_coverage.len()];
        return;
    }

    // This is the important part
    if bimbam {
        mw.matrix_f32 = vec![vec![0.0; len1]; pc.normalized_coverage.len()];
    } else {
        mw.matrix_bit =
            vec![BitVec::<u8, Lsb0>::repeat(false, len1 * 2); pc.normalized_coverage.len()];
    }
}

pub fn matrick_pack_wrapper(
    mw: &mut MatrixWrapper,
    pc: &mut PackCompact,
    want_node: bool,
    keep_zeros: bool,
    fraction: f32,
    method: Method,
    bimbam: bool,
    index_file: &Vec<u32>,
    index: usize,
    name: &String,
    absolute: u32,
) {
    // Add samples to matrix
    mw.sample_names.push(name.clone());

    if !index_file.is_empty() {
        pc.node_index = index_file.clone();
    }

    if pc.bin_coverage.is_empty() {
        if want_node {
            pc.calc_node_cov();
        } else if pc.normalized_coverage.is_empty() {
            pc.normalized_coverage = pc.coverage.iter().map(|x| *x as f32).collect();
        }

        let thresh;
        if absolute > 0 {
            thresh = absolute as f32;
        } else {
            thresh = PackCompact::get_threshold(pc, keep_zeros, fraction, 0.0, method);
        }

        if bimbam {
            f32_to_f32(mw, pc, thresh, index);
        } else {
            f32_to_bin(mw, pc, thresh, index);
        }
    } else {
        bin2bin(mw, pc, index);
    }
}

pub fn matrix_pack_wrapper(
    matrix_w: &mut MatrixWrapper,
    input: &Vec<PackCompact>,
    index: &Vec<u32>,
) {
    let first_entry = &input[0];

    if first_entry.data_type == DataType::TypeBit {
        matrix_w.matrix_bit =
            vec![BitVec::<u8, Lsb0>::repeat(false, input.len() * 2); first_entry.length as usize];

        for (i2, x) in input.iter().enumerate() {
            matrix_w.sample_names.push(x.name.clone());
            for (i, y) in x.bin_coverage.iter().enumerate() {
                if y == &true {
                    matrix_w.matrix_bit[i].get_mut(i2 * 2).unwrap().set(true);
                    matrix_w.matrix_bit[i]
                        .get_mut(i2 * 2 + 1)
                        .unwrap()
                        .set(true);
                }
            }
        }
        matrix_w.shape = (matrix_w.matrix_bit.len(), matrix_w.matrix_bit[0].len());
        matrix_w.sample_names = input.iter().map(|x| x.name.clone()).collect();
        if first_entry.is_sequence {
            matrix_w.geno_names = index2node_seq(index);
        } else {
            matrix_w.geno_names = remove_duplicates(index);
        }
    } else if first_entry.data_type == DataType::TypeU16 {
        let ll = first_entry.length as usize;
        matrix_w.matrix_u16 = vec![vec![0; input.len()]; ll];
        for (i, x) in input.iter().enumerate() {
            let d = &x.coverage;
            matrix_w.sample_names.push(x.name.clone());
            for (i2, y) in d.iter().enumerate() {
                matrix_w.matrix_u16[i2][i] = *y
            }
        }
        matrix_w.shape = (matrix_w.matrix_u16.len(), matrix_w.matrix_u16[0].len());
        matrix_w.sample_names = input.iter().map(|x| x.name.clone()).collect();

        if first_entry.is_sequence {
            matrix_w.geno_names = index2node_seq(index);
        } else {
            matrix_w.geno_names = remove_duplicates(index);
        }
    } else {
        let ll = first_entry.length as usize;
        matrix_w.matrix_f32 = vec![vec![0.0; input.len()]; ll];
        for (i, x) in input.iter().enumerate() {
            let d = &x.normalized_coverage;
            matrix_w.sample_names.push(x.name.clone());
            for (i2, y) in d.iter().enumerate() {
                matrix_w.matrix_f32[i2][i] = *y
            }
        }
        matrix_w.shape = (matrix_w.matrix_f32.len(), matrix_w.matrix_f32[0].len());
        matrix_w.sample_names = input.iter().map(|x| x.name.clone()).collect();

        if first_entry.is_sequence {
            matrix_w.geno_names = index2node_seq(index);
        } else {
            matrix_w.geno_names = remove_duplicates(index);
        }
    }
}

///
pub fn remove_duplicates(sorted_vec: &Vec<u32>) -> Vec<u64> {
    let mut unique_vec: Vec<u32> = sorted_vec.to_vec();
    let mut result = Vec::new();

    // If the input vector is empty, return an empty vector
    if sorted_vec.is_empty() {
        return result;
    }

    // Add the first element of the sorted vector
    unique_vec.dedup();

    for x in unique_vec.iter() {
        result.push(merge_u32_to_u64(*x, 0));
    }

    result
}
