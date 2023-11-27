use std::collections::{HashMap, HashSet};
use bitvec::order::Msb0;
use bitvec::vec::BitVec;
use gfa_reader::{GraphWrapper, NCGfa, NCPath};
use crate::core::core::MatrixWrapper;

pub fn gfa_nodes_reader(graph: &NCGfa<()>, gw: &Haplotype, wrapper: &GraphWrapper<NCPath>, matrix: &mut MatrixWrapper, bin: bool) {
    if bin {
        matrix.matrix_bin = vec![BitVec::<u8, Msb0>::repeat(false, gw.haplotype.len() * 2); graph.nodes.len()];
        matrix.shape = (matrix.matrix_bin.len(), matrix.matrix_bin[0].len());
        for (index, nn) in gw.haplotype.iter().enumerate() {
            println!("{}", nn.0);

            matrix.bin_names.push(nn.0.clone());
            matrix.names.push(nn.0.clone());
            for y in nn.1[0].1.iter(){
                for node in y.nodes.iter() {
                    matrix.matrix_bin[(node - 1) as usize].get_mut(index * 2).unwrap().set(true);
                }
            }
            for y in nn.1[1].1.iter(){
                for node in y.nodes.iter() {
                    matrix.matrix_bin[(node - 1) as usize].get_mut(index * 2 + 1).unwrap().set(true);
                }
            }

        }
    } else {
        matrix.matrix_core = vec![vec![0; wrapper.genomes.len()]; graph.nodes.len()];
        matrix.shape = (matrix.matrix_core.len(), matrix.matrix_core[0].len());
        for (index, node) in wrapper.genomes.iter().enumerate(){
            matrix.names.push(node.0.clone());

            for x in node.1.iter(){
                for node in x.nodes.iter() {
                    matrix.matrix_core[(node-1) as usize][index] += 1;
                }
            }
        }
    }
}



pub fn dir_nodes2(graph: &NCGfa<()>, gw: &Haplotype, wrapper: &GraphWrapper<NCPath>, matrix: &mut MatrixWrapper, bin: bool) {
    if bin {
        matrix.matrix_bin = vec![BitVec::<u8, Msb0>::repeat(false, gw.haplotype.len() * 2); graph.nodes.len()*2];
        matrix.shape = (matrix.matrix_bin.len(), matrix.matrix_bin[0].len());
        for (index, nn) in gw.haplotype.iter().enumerate() {

            matrix.names.push(nn.0.clone());
            for y in nn.1[0].1.iter(){
                for (node, dir) in y.nodes.iter().zip(y.dir.iter()) {
                    let i = (node-1) as usize * 2 + *dir as usize;
                    matrix.matrix_bin[i].get_mut(index * 2).unwrap().set(true);
                }
            }
            for y in nn.1[1].1.iter(){
                for (node, dir) in y.nodes.iter().zip(y.dir.iter()) {
                    let i = (node-1) as usize * 2 + *dir as usize;

                    matrix.matrix_bin[i].get_mut(index * 2 + 1).unwrap().set(true);
                }
            }

        }
    } else {
        matrix.matrix_core = vec![vec![0; wrapper.genomes.len()]; graph.nodes.len()*2];
        matrix.shape = (matrix.matrix_core.len(), matrix.matrix_core[0].len());
        for (index, node) in wrapper.genomes.iter().enumerate(){
            matrix.names.push(node.0.clone());

            for x in node.1.iter(){
                for (node, dir) in x.nodes.iter().zip(x.dir.iter()) {
                    let i = (node -1) as usize * 2 + *dir as usize;
                    matrix.matrix_core[i][index] += 1;
                }
            }
        }
    }
}

pub fn egdes1(graph: &NCGfa<()>, gw: &GraphWrapper<NCPath>, haplo: &Haplotype, matrix: &mut MatrixWrapper, bin: bool, names: &mut Vec<(u32, bool, u32, bool)>){
    let mut v: HashSet<(u32, bool, u32, bool)> = HashSet::new();
    for x in graph.edges.as_ref().unwrap().iter(){
        let j: u32 = x.to;
        let j2: u32 = x.from;
        v.insert((j2, x.from_dir, j, x.to_dir));
    }

    let mut k: Vec<_> = v.into_iter().collect();
    k.sort_by_key(|k| k.0);
    let mut bw = HashMap::new();
    for (index, x) in k.iter().enumerate(){

        bw.insert(x.clone(), index);
        names.push(x.clone());
    }
    if bin {
        matrix.matrix_bin = vec![BitVec::<u8, Msb0>::repeat(false, haplo.haplotype.len() * 2); bw.len()];
        matrix.shape = (matrix.matrix_bin.len(), matrix.matrix_bin[0].len());
        for (index, nn) in haplo.haplotype.iter().enumerate() {
            matrix.bin_names.push(nn.0.clone());
            matrix.names.push(nn.0.clone());
            for y in nn.1[0].1.iter() {
                for (o, node) in y.nodes.iter().zip(y.dir.iter()).enumerate().skip(1)  {
                    let ii = bw.get(&(y.nodes[o-1], y.dir[o-1], *node.0, *node.1)).unwrap();
                    matrix.matrix_bin[*ii].get_mut(index * 2).unwrap().set(true);
                }
            }
            for y in nn.1[1].1.iter() {
                for (o, node) in y.nodes.iter().zip(y.dir.iter()).enumerate().skip(1)  {
                    let ii = bw.get(&(y.nodes[o-1], y.dir[o-1], *node.0, *node.1)).unwrap();
                    matrix.matrix_bin[*ii].get_mut(index * 2 +1).unwrap().set(true);
                }
            }
        }
    } else {
        matrix.matrix_core = vec![vec![0; gw.genomes.len()]; bw.len()];
        matrix.shape = (matrix.matrix_core.len(), matrix.matrix_core[0].len());
        for (index, node) in gw.genomes.iter().enumerate(){
            matrix.names.push(node.0.clone());

            for y in node.1.iter(){
                for (o, node) in y.nodes.iter().zip(y.dir.iter()).enumerate().skip(1)  {
                    let ii = bw.get(&(y.nodes[o-1], y.dir[o-1], *node.0, *node.1)).unwrap();
                    *matrix.matrix_core[*ii].get_mut(index ).unwrap() += 1;
                }
            }
        }
    }
}







pub struct Haplotype<'a>{
    pub length: usize,
    pub haplotype: HashMap<String, [&'a (String, Vec<&'a NCPath>); 2]>,
    pub haplotype2: HashMap<String, [usize; 2]>,
}

impl <'a> Haplotype<'a> {
    pub fn from_wrapper(wrapper: &'a GraphWrapper<NCPath>, sep: &str) -> Self{
        let mut f: HashMap<String, [&'a (String, Vec<&'a NCPath>); 2]> = HashMap::new();
        let mut f2: HashMap<String, [usize; 2]> = HashMap::new();

        for (i, x) in wrapper.genomes.iter().enumerate(){
            let sstring: Vec<String> = x.0.split(sep).map(|s| s.to_string()).collect();
            let a = sstring[0].clone();
            if f.contains_key(&a){
                if f.get(&a).unwrap()[0].0 != f.get(&a).unwrap()[1].0{
                    panic!("daskdaj")
                } else {
                    f.get_mut(&a).unwrap()[1] = x;
                    f2.get_mut(&a).unwrap()[1] = i;

                }
            } else {
                f.insert(a.clone(), [x, x]);
                f2.insert(a, [i, i]);

            }

        }
        Self{
            length: 0,
            haplotype: f,
            haplotype2: f2,
        }

    }

    pub fn new() -> Self{
        Self{
            length: 0,
            haplotype: HashMap::new(),
            haplotype2: HashMap::new(),
        }
    }

    pub fn from_string_vec(vec1: Vec<String>, sep : &str) -> Self {
        let mut f2: HashMap<String, [usize; 2]> = HashMap::new();
        for (i, x) in vec1.iter().enumerate(){
            let sstring: Vec<String> = x.split(sep).map(|s| s.to_string()).collect();
            let a = sstring[0].clone();
            if f2.contains_key(&a){
                if f2.get(&a).unwrap()[0] != f2.get(&a).unwrap()[1]{
                    panic!("daskdaj")
                } else {
                    f2.get_mut(&a).unwrap()[1] = i;

                }
            } else {
                f2.insert(a, [i, i]);

            }

        }
        Self{
            length: 0,
            haplotype: HashMap::new(),
            haplotype2: f2,
        }
    }

}