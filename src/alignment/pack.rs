use crate::core::core::MatrixWrapper;
use crate::core::helper::index2node_seq;
use bitvec::order::Lsb0;
use bitvec::prelude::BitVec;
use packing_lib::core::core::{DataType, PackCompact};

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
    }
}

fn remove_duplicates(sorted_vec: &Vec<u32>) -> Vec<u64> {
    let mut unique_vec = Vec::new();

    // If the input vector is empty, return an empty vector
    if sorted_vec.is_empty() {
        return unique_vec;
    }

    // Add the first element of the sorted vector
    unique_vec.push(sorted_vec[0] as u64);

    // Iterate through the sorted vector and add only distinct elements
    for i in 1..sorted_vec.len() {
        // If the current element is different from the previous one, add it to the unique vector
        if sorted_vec[i] != sorted_vec[i - 1] {
            unique_vec.push(sorted_vec[i] as u64);
        }
    }

    unique_vec
}
