use crate::matrix::Matrix;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::io::{Write, BufWriter};
use std::fs::File;
use gfaR_wrapper::{GraphWrapper, NGfa};
use crate::helper::{binary2dec_bed, transpose, trans2};
use packing_lib::reader::{ReaderU16, wrapper_bool, ReaderBit, wrapper_u16, get_file_as_byte_vec};
use std::hash::Hash;
use bimap::{BiHashMap, BiMap};


/// Core data structure, which includes ever
pub struct MatrixWrapper<T: Debug>{
    pub matrix: Matrix,
    pub column_name: HashMap<u32, String>,
    pub row_name: HashMap<T, usize>,
    pub matrix_bin: Vec<Vec<bool>>,
    pub row2: BiMap<T, usize>,

}

impl <T>MatrixWrapper<T>

    where
        T: Debug + std::hash::Hash + std::cmp::Eq
{
    pub fn new() -> Self{
        let matrx = Matrix::new();
        let col: HashMap<u32, String> = HashMap::new();
        let row: HashMap<T, usize> = HashMap::new();
        let matrx2 = Vec::new();
        let row2 = BiMap::new();
        Self{
            matrix: matrx,
            column_name: col,
            row_name: row,
            matrix_bin: matrx2,
            row2: row2,
        }
    }

    /// Writing bim file
    /// Information: https://www.cog-genomics.org/plink/1.9/formats#bim
    pub fn write_bim(&self, out_prefix: &str, t: &str){
        let f = File::create([out_prefix, t,  "bim"].join(".")).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        let mut helper_vec = Vec::new();
        for (k,v) in self.row_name.iter(){
            helper_vec.push((v,k));
        }
        helper_vec.sort_by_key(|k| k.0);



        for (k,_v) in helper_vec.iter(){
            write!(f, "{}\t{}\t{}\t{}\t{}\t{}\n", "graph", ".", 0, k, "A", "T").expect("Not able to write ");
        }
        let f = File::create([out_prefix, t,  "bimhelper"].join(".")).expect("Unable to create file");
        let mut f = BufWriter::new(f);
        for (v,k) in helper_vec.iter(){
            write!(f, "{}\t{:?}\n", v, k).expect("Not able to write");

        }
    }

    /// Writing BED File
    pub fn write_bed(&self, out_prefix: &str, t: &str){
        //hexdump -C test.bin
        // xxd -b file
        // xxd file
        // SNP: 00000001 , 0
        // IND: 00000000, 1

        println!("{} {}", self.matrix_bin.len(), self.matrix_bin[0].len());
        let mut buff: Vec<u8> = vec![108, 27, 1];
        let h2 = trans2( &self.matrix_bin);
        println!("{} {}", h2.len(), h2[0].len());
        for x in h2.iter(){
            let j: Vec<&[&bool]> = x.chunks(4).collect();
            for x in j{
                buff.push(binary2dec_bed(x));
            }
            //println!("Number of bytes {}", buff.len());
            //println!("x {}", x.len());
        }


        let mut file = File::create([out_prefix, t, "bed"].join(".")).expect("Not able to write ");
        file.write_all(&buff).expect("Not able to write ")

    }

    /// Writing file wrapper
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
            write!(f, "{}\n", self.column_name.get(&(x as u32)).unwrap()).expect("Can not write file");
        }
    }

    pub fn remove_stuff(& mut self, v: Vec<usize>) {
        for x in v.iter(){
            self.row2.remove_by_right(x);
        }
    }



}

//-------------------------------------------
// Modification

pub fn mod_row(mu: & mut MatrixWrapper<u32>, v: Vec<usize>){
    for x in v.iter(){
        mu.row2.remove_by_right(&1);
    }

}

pub fn mod_row2(mu: & mut MatrixWrapper<(u32, bool)>, v: Vec<usize>){
    for x in v.iter(){
        mu.row2.remove_by_right(&1);
    }

}

pub fn mod_row3(mu: & mut MatrixWrapper<(u32, bool, u32, bool)>, v: Vec<usize>){
    for x in v.iter(){
        mu.row2.remove_by_right(&1);
    }

}



//--------------------------------------------
// Construction of Matrix-Wrapper

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
        eprintln!("{}", name);
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

    let mut v: HashSet<(u32, bool)> = HashSet::new();

    for x in graph.paths.iter(){
        for x2 in 0..x.nodes.len(){
            v.insert((x.nodes[x2], x.dir[x2]));
        }
    }
    let mut k: Vec<_> = v.into_iter().collect();
    k.sort_by_key(|k| k.0);
    let mut count = 0;
    for x in k.iter(){
        mw.row_name.insert(x.clone(), count);
        count += 1;
    }
    let mut bim: BiHashMap<(u32, bool), usize> = bimap::BiMap::new();
    let mut count = 0;
    for x in k.iter(){
        bim.insert(x.clone(), count);
        mw.row_name.insert(x.clone(), count);
        count += 1;
    }
    //eprintln!("{:?}", bim);







    for (index, (name, vec1)) in gwrapper.genomes.iter().enumerate(){
        mw.column_name.insert(index as u32, name.clone());
        //println!("{}", name);
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
    let mut v: HashSet<(u32, bool, u32, bool)> = HashSet::new();
    for (i, x) in graph.edges.iter().enumerate(){
        let j: u32 = x.to;
        let j2: u32 = x.from;
        v.insert((j2, x.from_dir, j, x.to_dir));
    }

    let mut k: Vec<_> = v.into_iter().collect();
    k.sort_by_key(|k| k.0);
    let mut count = 0;
    for x in k.iter(){
        mw.row_name.insert(x.clone(), count);
        count += 1;
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


/// Make matrix from bit vector
pub fn matrix_node_coverage(filename: &str) -> MatrixWrapper<u32>{
    let g: Vec<u8> = get_file_as_byte_vec(filename);
    let k: Vec<ReaderBit> = wrapper_bool(&g);

    let mut mw: MatrixWrapper<u32>  = MatrixWrapper::new();
    for (i,x) in k.iter().enumerate(){
        //println!("{}", x.name);
        mw.column_name.insert(i as u32,x.name.clone());
        let mut k : Vec<bool> = Vec::new();
        for y in x.cc.iter(){

            k.push(*y);
        }
        //println!("{}", k.len());
        mw.matrix_bin.push(k);
    }
    for x in 0..mw.matrix_bin[0].len(){
        mw.row_name.insert(x as u32, x);
    }
    mw
}


pub fn matrix_node_coverage2(filename: &str) -> MatrixWrapper<u32>{
    let g: Vec<u8> = get_file_as_byte_vec(filename);
    let k: Vec<ReaderU16> = wrapper_u16(&g);

    let mut mw: MatrixWrapper<u32>  = MatrixWrapper::new();
    for (i,x) in k.iter().enumerate(){
        //println!("{}", x.name);
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
