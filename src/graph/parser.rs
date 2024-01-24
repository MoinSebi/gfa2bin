use crate::core::core::MatrixWrapper;
use crate::core::helper::{merge_u32_to_u64, to_string1, Feature};
use bitvec::macros::internal::funty::Fundamental;
use bitvec::order::Lsb0;
use bitvec::vec::BitVec;
use gfa_reader::{NCPath, Pansn};

pub fn gfa_reader(
    matrix: &mut MatrixWrapper,
    graph_wrapper: &Pansn<NCPath>,
    bin: bool,
    feature: Feature,
) {
    if bin {
        matrix.matrix_bin = vec![
            BitVec::<u8, Lsb0>::repeat(false, graph_wrapper.genomes.len() * 2);
            matrix.geno_names.len()
        ];
        matrix.shape = (matrix.matrix_bin.len(), matrix.matrix_bin[0].len());
        let hm = &matrix.geno_names;
        for (path_index, nn) in graph_wrapper.genomes.iter().enumerate() {
            matrix.sample_names.push(nn.name.clone());
            if nn.haplotypes.len() == 1 {
                for ncpath in nn.haplotypes[0].paths.iter() {
                    let mut iter1 = paths_to_u64vec(ncpath, feature);
                    iter1.sort();
                    let mut i = 0;
                    let mut j = 0;

                    while i < iter1.len() && j < hm.len() {
                        if iter1[i] == hm[j] {
                            matrix.matrix_bin[j]
                                .get_mut(path_index * 2 + 1)
                                .unwrap()
                                .set(true);

                            matrix.matrix_bin[j]
                                .get_mut(path_index * 2)
                                .unwrap()
                                .set(true);

                            i += 1;
                            j += 1;
                        } else if iter1[i] < hm[j] {
                            i += 1;
                        } else {
                            j += 1;
                        }
                    }
                }
            } else if nn.haplotypes.len() == 2 {
                for ncpath in nn.haplotypes[0].paths.iter() {
                    let mut iter1 = paths_to_u64vec(ncpath, feature);
                    iter1.sort();
                    let mut i = 0;
                    let mut j = 0;

                    while i < iter1.len() && j < hm.len() {
                        if iter1[i] == hm[j] {
                            matrix.matrix_bin[j]
                                .get_mut(path_index * 2)
                                .unwrap()
                                .set(true);
                            i += 1;
                            j += 1;
                        } else if iter1[i] < hm[j] {
                            i += 1;
                        } else {
                            j += 1;
                        }
                    }
                }
                for ncpath in nn.haplotypes[1].paths.iter() {
                    let iter1 = paths_to_u64vec(ncpath, feature);
                    let mut i = 0;
                    let mut j = 0;

                    while i < iter1.len() && j < hm.len() {
                        if iter1[i] == hm[j] {
                            let check = matrix.matrix_bin[j][path_index * 2];
                            if check {
                                matrix.matrix_bin[j]
                                    .get_mut(path_index * 2 + 1)
                                    .unwrap()
                                    .set(true);
                            } else {
                                matrix.matrix_bin[j]
                                    .get_mut(path_index * 2)
                                    .unwrap()
                                    .set(true);
                            }
                            i += 1;
                            j += 1;
                        } else if iter1[i] < hm[j] {
                            i += 1;
                        } else {
                            j += 1;
                        }
                    }
                }
            } else {
                panic!("Not implemented");
            }
        }
    }
}

pub fn paths_to_u64vec(path: &NCPath, feature: Feature) -> Vec<u64> {
    let mut vec_u64 = Vec::new();
    for i in 0..path.nodes.len() - 1 {
        let v1 = path.nodes[i];
        let v2 = path.dir[i];
        let v3 = path.nodes[i + 1];
        let v4 = path.dir[i + 1];

        if feature == Feature::Node {
            vec_u64.push(v1 as u64);
        } else if feature == Feature::DirNode {
            vec_u64.push(v1 as u64 * 2 + v2 as u64);
        } else if feature == Feature::Edge {
            let u1 = v1 * 2 + v2 as u32;
            let u2 = v3 * 2 + v4 as u32;
            vec_u64.push(merge_u32_to_u64(u1, u2));
        }
    }
    if feature == Feature::Node {
        vec_u64.push(path.nodes[path.nodes.len() - 1] as u64);
    } else if feature == Feature::DirNode {
        vec_u64.push(
            path.nodes[path.nodes.len() - 1] as u64 * 2 + path.dir[path.nodes.len() - 1] as u64,
        );
    }
    vec_u64
}
