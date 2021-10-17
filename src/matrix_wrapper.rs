use crate::matrix::Matrix;
use std::collections::HashMap;
use std::fmt::Debug;
use std::io::{Write, BufWriter};
use std::fs::File;
use gfaR_wrapper::{GraphWrapper, NGfa};
use crate::pack_reader::{read_in2, read_in3};
use crate::helper::{binary2dec_bed, transpose};

pub struct MatrixWrapper<T: Debug>{
    pub matrix: Matrix,
    pub column_name: HashMap<u32, String>,
    pub row_name: HashMap<T, usize>,
    pub matrix_bin: Vec<Vec<bool>>,

}

impl <T>MatrixWrapper<T>

    where
        T: Debug
{
    pub fn new() -> Self{
        let matrx = Matrix::new();
        let col: HashMap<u32, String> = HashMap::new();
        let row: HashMap<T, usize> = HashMap::new();
        let matrx2 = Vec::new();
        Self{
            matrix: matrx,
            column_name: col,
            row_name: row,
            matrix_bin: matrx2,
        }
    }

    /// Wirting bim file
    /// Information: https://www.cog-genomics.org/plink/1.9/formats#bim
    pub fn write_bim(&self, out_prefix: &str, t: &str){
        let f = File::create([out_prefix, t,  "bim"].join(".")).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        let mut helper_vec = Vec::new();
        for (k,v) in self.row_name.iter(){
            helper_vec.push((v,k));
        }
        helper_vec.sort_by_key(|k| k.0);



        for (k,v) in helper_vec.iter(){
            write!(f, "{}\t{}\t{}\t{}\t{}\t{}\n", "graph", ".", 0, k, "A", "T").expect("Not able to write ");
        }
        let f = File::create([out_prefix, t,  "bimhelper"].join(".")).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        for (v,k) in helper_vec.iter(){
            write!(f, "{}\t{:?}\n", v, k).expect("Not able to write");

        }
    }

    pub fn write_bed(&self, out_prefix: &str, t: &str){
        //hexdump -C test.bin
        // xxd -b file
        // xxd file
        // SNP: 00000001 , 0
        // IND: 00000000, 1

        println!("{} {}", self.matrix_bin.len(), self.matrix_bin[0].len());
        let mut buff: Vec<u8> = vec![108, 27, 1];
        let h2 = transpose( &self.matrix_bin);
        println!("{} {}", h2.len(), h2[0].len());
        for x in h2.iter(){
            let j: Vec<&[bool]> = x.chunks(4).collect();
            for x in j{
                buff.push(binary2dec_bed(x));
            }
            //println!("Number of bytes {}", buff.len());
            //println!("x {}", x.len());
        }


        let mut file = File::create([out_prefix, t, "bed"].join(".")).expect("Not able to write ");
        file.write_all(&buff).expect("Not able to write ")

    }

    pub fn write(& mut self, what: &str, out_prefix: &str, t: &str){
        if (what == "bed") | (what == "all"){
            if self.matrix_bin.is_empty(){
                self.matrix_bin = self.matrix.copy(1);
            };
            self.write_bed(out_prefix, t);
            self.write_bim(out_prefix, t);
            self.write_helper2(out_prefix);
        }
        if (what == "bimbam") | (what == "all"){
            self.matrix.write_bimbam(out_prefix, t);
        }
    }

    pub fn write_helper2(&self, out_prefix: &str){
        let f = File::create([out_prefix,  "bimbim"].join(".")).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        for x in 0..self.column_name.len(){
            write!(f, "{}\n", self.column_name.get(&(x as u32)).unwrap());
        }
    }

}


/// Make matrix for nodes
pub fn matrix_node(gwrapper: &GraphWrapper, graph: &NGfa) -> MatrixWrapper<u32>{
    let mut mw: MatrixWrapper<u32> = MatrixWrapper::new();
    let mut h: Vec<u32> = graph.nodes.keys().cloned().collect();

    h.sort_by(|a, b| a.partial_cmp(b).unwrap());

    for (i, x) in h.iter().enumerate(){
        let node:u32 = x.clone();
        mw.row_name.insert (node, i);
    }


    for (index, (name, paths)) in gwrapper.genomes.iter().enumerate(){
        println!("{}", name);
        mw.column_name.insert( index as u32, name.clone());
        let mut nody: Vec<u32> = vec![0; mw.row_name.len()] ;
        for x in paths.iter(){
            for y in x.nodes.iter(){
                nody[mw.row_name[&y]] += 1;

            }
        }
        mw.matrix.matrix_core.push(nody);
    }
    mw

}

/// Make matrix for directed nodes
pub fn matrix_dir_node(gwrapper: &GraphWrapper, graph: &NGfa) -> MatrixWrapper<(u32, bool)>{
    let mut mw: MatrixWrapper<(u32, bool)> = MatrixWrapper::new();
    let mut count = 0;
    for x in graph.paths.iter(){
        for x2 in 0..x.nodes.len(){
            let node: u32 = x.nodes[x2].clone();
            if ! mw.row_name.contains_key(&(node, x.dir[x2])){
                mw.row_name.insert((node, x.dir[x2].clone()), count);
                count += 1;
            }
        }
    }
    for (index, (name, vec1)) in gwrapper.genomes.iter().enumerate(){
        mw.column_name.insert(index as u32, name.clone());
        println!("{}", name);
        let mut dir_nodes : Vec<u32> = vec![0; mw.row_name.len()] ;
        for x in vec1.iter(){
            for x2 in 0..x.nodes.len(){
                let node: u32 = x.nodes[x2].clone();

                dir_nodes[mw.row_name[&(node, x.dir[x2].clone())]] += 1;

            }
        }
        mw.matrix.matrix_core.push(dir_nodes);
    }
    mw
}

/// Make matrix for edges
pub fn matrix_edge(gwrapper: &GraphWrapper, graph: &NGfa) -> MatrixWrapper<(u32, bool, u32, bool)>{
    let mut mw: MatrixWrapper<(u32, bool, u32, bool)> = MatrixWrapper::new();
    for (i, x) in graph.edges.iter().enumerate(){
        let j: u32 = x.to;
        let j2: u32 = x.from;
        mw.row_name.insert((j2, x.from_dir, j, x.to_dir), i);
    }
    for (index, (name, vec1)) in gwrapper.genomes.iter().enumerate(){
        mw.column_name.insert( index as u32, name.clone());
        let mut dir_nodes : Vec<u32> = vec![0; mw.row_name.len()] ;
        for x in vec1.iter(){
            for x2 in 0..x.nodes.len()-1{
                let node: u32 = x.nodes[x2];
                let node2: u32 = x.nodes[x2+1];
                dir_nodes[mw.row_name[&(node, x.dir[x2], node2, x.dir[x2+1])]] += 1;

            }
        }
        mw.matrix.matrix_core.push(dir_nodes);
    }
    mw
}


/// Make matrix for edges
pub fn matrix_node_coverage(cov: Vec<read_in2>) -> MatrixWrapper<u32>{
    let mut mw: MatrixWrapper<u32>  = MatrixWrapper::new();
    for (i,x) in cov.iter().enumerate(){
        println!("{}", x.name);
        mw.column_name.insert(i as u32,x.name.clone());
        let mut k : Vec<bool> = Vec::new();
        for y in x.cc.iter(){

            k.push(*y);
        }
        println!("{}", k.len());
        mw.matrix_bin.push(k);
    }
    for x in 0..mw.matrix_bin[0].len(){
        mw.row_name.insert(x as u32, x);
    }
    mw
}


pub fn matrix_node_coverage2(cov: Vec<read_in3>) -> MatrixWrapper<u32>{
    let mut mw: MatrixWrapper<u32>  = MatrixWrapper::new();
    for (i,x) in cov.iter().enumerate(){
        println!("{}", x.name);
        mw.column_name.insert(i as u32,x.name.clone());
        let mut k : Vec<u32> = Vec::new();
        for y in x.cc.iter(){

            k.push(*y as u32);
        }
        mw.matrix.matrix_core.push(k);
    }
    for x in 0..mw.matrix.matrix_core[0].len(){
        mw.row_name.insert(x as u32, x);
    }
    mw
}
