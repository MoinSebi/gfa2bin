use std::fs::File;
use std::io::{BufRead, BufReader};
use crate::core::core::MatrixWrapper;
use crate::core::helper::{merge_u32_to_u64, Feature};

use bitvec::order::Lsb0;
use bitvec::vec::BitVec;
use gfa_reader::{Pansn, Path};


/// Read a gfa file and convert it to a matrix (bit or u16)
///
/// Read the graph sample by
pub fn gfa_reader(
    matrix: &mut MatrixWrapper,
    graph_wrapper: &Pansn<u32, (), ()>,
    want_bool: bool,
    feature: Feature,
) {
    // Don't count
    if want_bool {
        // Matrix bit
        matrix.matrix_bit = vec![
            BitVec::<u8, Lsb0>::repeat(false, graph_wrapper.genomes.len() * 2);
            matrix.geno_names.len()
        ];
        matrix.shape = (matrix.matrix_bit.len(), matrix.matrix_bit[0].len());

        let index2geno = &matrix.geno_names;
        for (path_index, nn) in graph_wrapper.genomes.iter().enumerate() {
            matrix.sample_names.push(nn.name.clone());

            // "haploid"
            if nn.haplotypes.len() == 1 {
                for ncpath in nn.haplotypes[0].paths.iter() {
                    let path_geno_vec = paths_to_u64vec(ncpath, feature);
                    let mut i = 0;
                    let mut j = 0;

                    while i < path_geno_vec.len() && j < index2geno.len() {
                        if path_geno_vec[i] == index2geno[j] {
                            matrix.matrix_bit[j]
                                .get_mut(path_index * 2 + 1)
                                .unwrap()
                                .set(true);

                            matrix.matrix_bit[j]
                                .get_mut(path_index * 2)
                                .unwrap()
                                .set(true);

                            i += 1;
                            j += 1;
                        } else if path_geno_vec[i] < index2geno[j] {
                            i += 1;
                        } else {
                            j += 1;
                        }
                    }
                }

                // "diploid"
            } else if nn.haplotypes.len() == 2 {
                // First haplotypes
                for ncpath in nn.haplotypes[0].paths.iter() {
                    let iter1 = paths_to_u64vec(ncpath, feature);
                    let mut i = 0;
                    let mut j = 0;

                    while i < iter1.len() && j < index2geno.len() {
                        if iter1[i] == index2geno[j] {
                            matrix.matrix_bit[j]
                                .get_mut(path_index * 2)
                                .unwrap()
                                .set(true);
                            i += 1;
                            j += 1;
                        } else if iter1[i] < index2geno[j] {
                            i += 1;
                        } else {
                            j += 1;
                        }
                    }
                }

                // Second haplotypes
                for ncpath in nn.haplotypes[1].paths.iter() {
                    let iter1 = paths_to_u64vec(ncpath, feature);
                    let mut i = 0;
                    let mut j = 0;

                    while i < iter1.len() && j < index2geno.len() {
                        if iter1[i] == index2geno[j] {
                            let check = matrix.matrix_bit[j][path_index * 2];
                            if check {
                                matrix.matrix_bit[j]
                                    .get_mut(path_index * 2 + 1)
                                    .unwrap()
                                    .set(true);
                            } else {
                                matrix.matrix_bit[j]
                                    .get_mut(path_index * 2)
                                    .unwrap()
                                    .set(true);
                            }
                            i += 1;
                            j += 1;
                        } else if iter1[i] < index2geno[j] {
                            i += 1;
                        } else {
                            j += 1;
                        }
                    }
                }
                // If you have more than 3 haplotypes
            } else {
                panic!("Not implemented");
            }
        }
        // Now count
    } else {
        matrix.matrix_u16 =
            vec![vec![0; graph_wrapper.get_haplo_path().len()]; matrix.geno_names.len()];
        matrix.shape = (matrix.matrix_u16.len(), matrix.matrix_u16[0].len());
        let hm = &matrix.geno_names;
        let mut c = 0;
        for nn in graph_wrapper.genomes.iter() {
            matrix.sample_names.push(nn.name.clone());
            if nn.haplotypes.len() == 1 {
                matrix.sample_index_u16.push([c, c]);
            } else if nn.haplotypes.len() == 2 {
                matrix.sample_index_u16.push([c, c + 1]);
            } else {
                panic!("Not implemented");
            }

            for haplotype in nn.haplotypes.iter() {
                for path in haplotype.paths.iter() {
                    let iter1 = paths_to_u64vec(path, feature);
                    let mut i = 0;
                    let mut j = 0;

                    let mut lastentry = 0;

                    while i < iter1.len() && j < hm.len() {
                        if iter1[i] == hm[j] {
                            matrix.matrix_u16[j][c] += 1;
                            lastentry = iter1[i];
                            i += 1;
                            j += 1;
                        } else if iter1[i] < hm[j] {
                            if iter1[i] == lastentry {
                                matrix.matrix_u16[j][c] += 1;
                            }
                            lastentry = iter1[i];
                            i += 1;
                        } else {
                            j += 1;
                        }
                    }
                }
                c += 1
            }
        }
    }
}

/// Convert a path to a vector of u64 depending on the feature
///
/// Node: Node ID
/// DirNode: Node ID + Direction
/// Edge: (Node1 ID, Direction 1) + (Node2 ID, Direction 2)
pub fn paths_to_u64vec(path: &Path<u32, (), ()>, feature: Feature) -> Vec<u64> {
    let mut vec_u64 = Vec::new();
    for i in 0..path.nodes.len() - 1 {
        let n1 = path.nodes[i];
        let d1 = path.dir[i];
        let n2 = path.nodes[i + 1];
        let d2 = path.dir[i + 1];

        if feature == Feature::Node {
            vec_u64.push(n1 as u64);
        } else if feature == Feature::DirNode {
            vec_u64.push(n1 as u64 * 2 + d1 as u64);
        } else if feature == Feature::Edge {
            let u1 = n1 * 2 + d1 as u32;
            let u2 = n2 * 2 + d2 as u32;
            // Merge u32 to u64
            vec_u64.push(merge_u32_to_u64(u1, u2));
        }
    }
    // Add the last entry
    if feature == Feature::Node {
        vec_u64.push(path.nodes[path.nodes.len() - 1] as u64);
    } else if feature == Feature::DirNode {
        vec_u64.push(
            path.nodes[path.nodes.len() - 1] as u64 * 2 + path.dir[path.nodes.len() - 1] as u64,
        );
    }
    vec_u64.sort();
    vec_u64
}
