use crate::core::core::MatrixWrapper;
use bitvec::order::Lsb0;
use bitvec::prelude::BitVec;
use packing_lib::core::core::PackCompact;

pub fn matrix_pack_wrapper(matrix_w: &mut MatrixWrapper, input: &Vec<PackCompact>, index: &Vec<u32>) {
    let tt = &input[0];

    if tt.is_binary {
        matrix_w.matrix_bin =
            vec![BitVec::<u8, Lsb0>::repeat(false, input.len() * 2); tt.length as usize];
        for (i2, x) in input.iter().enumerate() {
            matrix_w.sample_names.push(x.name.clone());
            for (i, y) in x.bin_coverage.iter().enumerate() {
                if y == &true {
                    matrix_w.matrix_bin[i].get_mut(i2 * 2).unwrap().set(true);
                    matrix_w.matrix_bin[i]
                        .get_mut(i2 * 2 + 1)
                        .unwrap()
                        .set(true);
                }
            }
        }
        matrix_w.shape = (matrix_w.matrix_bin.len(), matrix_w.matrix_bin[0].len());
        matrix_w.sample_names = input.iter().map(|x| x.name.clone()).collect();
        if tt.is_sequence {
            matrix_w.geno_names = index.into_iter().map(|n| *n as u64).collect();
        } else {
            let i2 = remove_duplicates(index);
            matrix_w.geno_names = i2.into_iter().map(|n| n as u64).collect();
        }

    } else {
        let ll = tt.length as usize;
        matrix_w.matrix_core = vec![vec![0; input.len()]; ll];
        for (i, x) in input.iter().enumerate() {
            let mut d = &vec![];

            if x.is_sequence {
                d = &x.coverage;
            } else {
                d = &x.node_coverage;
            }
            matrix_w.sample_names.push(x.name.clone());
            for (i2, y) in x.coverage.iter().enumerate() {
                matrix_w.matrix_core[i2][i] = *y as u32;
            }
        }
        matrix_w.shape = (matrix_w.matrix_core.len(), matrix_w.matrix_core[0].len());
        matrix_w.sample_names = input.iter().map(|x| x.name.clone()).collect();
        if tt.is_sequence {
            matrix_w.geno_names = index.into_iter().map(|n| *n as u64).collect();
        } else {
            let i2 = remove_duplicates(index);
            matrix_w.geno_names = i2.into_iter().map(|n| n as u64).collect();
        }
    }
}

fn remove_duplicates(sorted_vec: &Vec<u32>) -> Vec<u32> {
    let mut unique_vec = Vec::new();

    // If the input vector is empty, return an empty vector
    if sorted_vec.is_empty() {
        return unique_vec;
    }

    // Add the first element of the sorted vector
    unique_vec.push(sorted_vec[0]);

    // Iterate through the sorted vector and add only distinct elements
    for i in 1..sorted_vec.len() {
        // If the current element is different from the previous one, add it to the unique vector
        if sorted_vec[i] != sorted_vec[i - 1] {
            unique_vec.push(sorted_vec[i]);
        }
    }

    unique_vec
}