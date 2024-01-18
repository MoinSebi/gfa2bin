use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use bitvec::macros::internal::funty::Fundamental;
use bitvec::order::{Lsb0, Msb0};
use bitvec::vec::BitVec;
use gfa_reader::{Pansn, NCGfa, NCPath};
use crate::core::core::{MatrixWrapper};
use crate::core::helper::{Feature, GenoName};

pub fn gfa_nodes_reader2(matrix: &mut MatrixWrapper, graph_wrapper: &Pansn<NCPath>, graph: &NCGfa<()>, bin: bool, feature: &Feature){
    if bin {
        matrix.matrix_bin = vec![BitVec::<u8, Lsb0>::repeat(false, graph_wrapper.genomes.len() * 2); graph.nodes.len()];
        matrix.shape = (matrix.matrix_bin.len(), matrix.matrix_bin[0].len());
        let mut hm = &matrix.geno_map;
        for (path_index, nn) in graph_wrapper.genomes.iter().enumerate() {
            matrix.sample_names.push(nn.name.clone());
            if nn.haplotypes.len() == 1{
                for haplotype in nn.haplotypes.iter() {
                    for node1 in haplotype.paths.iter() {
                        for node in node1.nodes.iter() {
                            let index = hm.get(&GenoName{name: *node as u64}).unwrap();
                            matrix.matrix_bin[*index].get_mut(path_index * 2 + 1).unwrap().set(true);
                            matrix.matrix_bin[*index].get_mut(path_index * 2).unwrap().set(true);

                        }
                    }
                }
            } else if nn.haplotypes.len() == 2{
                for node1 in nn.haplotypes[0].paths.iter() {
                    for node in node1.nodes.iter(){
                        let index = hm.get(&GenoName{name: *node as u64}).unwrap();
                        matrix.matrix_bin[*index].get_mut(path_index * 2 ).unwrap().set(true);
                        }
                }
                for node1 in nn.haplotypes[1].paths.iter() {
                    for node in node1.nodes.iter(){
                        let index = hm.get(&GenoName{name:*node as u64}).unwrap();
                        let check = matrix.matrix_bin[*index][path_index*2];
                        if check.as_bool(){
                            matrix.matrix_bin[*index].get_mut(path_index * 2 +1).unwrap().set(true);
                        } else {
                            matrix.matrix_bin[*index].get_mut(path_index * 2 ).unwrap().set(true);
                        }
                    }
                }

            } else {
                panic!("Not implemented");
            }
        }
    }
}



