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
            matrix.geno_map.len()
        ];
        matrix.shape = (matrix.matrix_bin.len(), matrix.matrix_bin[0].len());
        let hm = &matrix.geno_map;
        for (path_index, nn) in graph_wrapper.genomes.iter().enumerate() {
            matrix.sample_names.push(nn.name.clone());
            if nn.haplotypes.len() == 1 {
                for ncpath in nn.haplotypes[0].paths.iter() {
                    let iter1 = paths_to_u64vec(ncpath, feature);
                    for node in iter1.iter() {
                        println!("{}", to_string1(*node, &feature));
                        let index = hm.get(node).unwrap();
                        matrix.matrix_bin[*index]
                            .get_mut(path_index * 2 + 1)
                            .unwrap()
                            .set(true);
                        matrix.matrix_bin[*index]
                            .get_mut(path_index * 2)
                            .unwrap()
                            .set(true);
                    }
                }
            } else if nn.haplotypes.len() == 2 {
                for node1 in nn.haplotypes[0].paths.iter() {
                    let iter1 = paths_to_u64vec(node1, feature);
                    for node in iter1.iter() {
                        let index = hm.get(node).unwrap();
                        matrix.matrix_bin[*index]
                            .get_mut(path_index * 2)
                            .unwrap()
                            .set(true);
                    }
                }
                for node1 in nn.haplotypes[1].paths.iter() {
                    let iter1 = paths_to_u64vec(node1, feature);
                    for node in iter1.iter() {
                        let index = hm.get(node).unwrap();
                        let check = matrix.matrix_bin[*index][path_index * 2];
                        if check.as_bool() {
                            matrix.matrix_bin[*index]
                                .get_mut(path_index * 2 + 1)
                                .unwrap()
                                .set(true);
                        } else {
                            matrix.matrix_bin[*index]
                                .get_mut(path_index * 2)
                                .unwrap()
                                .set(true);
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
    let mut v = Vec::new();
    for i in 0..path.nodes.len() - 1 {
        let v1 = path.nodes[i];
        let v2 = path.dir[i];
        let v3 = path.nodes[i + 1];
        let v4 = path.dir[i + 1];

        if feature == Feature::Node {
            v.push(v1 as u64);
        } else if feature == Feature::DirNode {
            v.push(v1 as u64 * 2 + v2 as u64);
        } else if feature == Feature::Edge {
            let u1 = v1 * 2 + v2 as u32;
            let u2 = v3 * 2 + v4 as u32;
            v.push(merge_u32_to_u64(u1, u2));
        }
    }
    if feature == Feature::Node {
        v.push(path.nodes[path.nodes.len() - 1] as u64);
    } else if feature == Feature::DirNode {
        v.push(path.nodes[path.nodes.len() - 1] as u64 * 2 + path.dir[path.nodes.len() - 1] as u64);
    }
    v
}
