use crate::matrix::Matrix;
use std::collections::HashMap;
use std::fmt::Debug;
use std::io::{Write, BufWriter};
use std::fs::File;
use gfaR_wrapper::{GraphWrapper, NGfa};

pub struct MatrixWrapper<T: Debug>{
    pub matrix: Matrix,
    pub column_name: HashMap<String, u32>,
    pub row_name: HashMap<T, usize>,

}

impl <T>MatrixWrapper<T>

    where
        T: Debug
{
    pub fn new() -> Self{
        let matrx = Matrix::new();
        let col: HashMap<String, u32> = HashMap::new();
        let row: HashMap<T, usize> = HashMap::new();
        Self{
            matrix: matrx,
            column_name: col,
            row_name: row,
        }
    }

    /// Wirting bim file
    /// Information: https://www.cog-genomics.org/plink/1.9/formats#bim
    pub fn write_bim(&self, out_prefix: &str, t: &str){
        let f = File::create([out_prefix, t,  "bim"].join(".")).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        let mut helper_vec = Vec::new();
        for (index, (k,v)) in self.row_name.iter().enumerate(){
            write!(f, "{}\t{}\t{}\t{}\t{}\t{}\n", "graph", ".", 0, index, "A", "T").expect("Not able to write ");
            helper_vec.push(k);
        }
        let f = File::create([out_prefix, t,  "bimhelper"].join(".")).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        for (i, x) in helper_vec.iter().enumerate(){
            write!(f, "{}\t{:?}\n", i, x).expect("Not able to write");

        }
    }

    pub fn write(&self, what: &str, out_prefix: &str, t: &str){
        if (what == "bed") | (what == "all"){
            self.matrix.write_bed(out_prefix, t);
            self.write_bim(out_prefix, t);
        }
        if (what == "bimbam") | (what == "all"){
            self.matrix.write_bimbam(out_prefix, t);
        }
    }

}


/// Make matrix for nodes
pub fn matrix_node(gwrapper: &GraphWrapper, graph: &NGfa) -> MatrixWrapper<u32>{
    let mut mw: MatrixWrapper<u32> = MatrixWrapper::new();


    for (i, x) in graph.nodes.iter().enumerate(){
        let node:u32 = x.0.clone();
        mw.row_name.insert (node, i);
    }


    for (index, (name, paths)) in gwrapper.genomes.iter().enumerate(){

        mw.column_name.insert(name.clone(), index as u32);
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
        mw.column_name.insert(name.clone(), index as u32);
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
        mw.column_name.insert(name.clone(), index as u32);
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