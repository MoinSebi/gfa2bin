use crate::matrix::Matrix;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::io::{Write, BufWriter, BufReader, BufRead};
use std::fs::File;
use gfaR_wrapper::{GraphWrapper, NGfa};
use crate::helper::{binary2dec_bed, trans2};
use packing_lib::reader::{ReaderU16, wrapper_bool, ReaderBit, wrapper_u16, get_file_as_byte_vec};
use bimap::{BiHashMap, BiMap};


/// Core data structure, which includes ever
pub struct MatrixWrapper2{
    pub matrix: Matrix,
    pub column_name: HashMap<u32, String>,
    pub matrix_bin: Vec<Vec<bool>>

}

impl MatrixWrapper2{

    /// Dummy initialization
    pub fn new() -> Self {
        let matrx = Matrix::new();
        let col: HashMap<u32, String> = HashMap::new();
        let matrx2: Vec<Vec<bool>> = Vec::new();
        Self {
            matrix: matrx,
            column_name: col,
            matrix_bin: matrx2,
        }
    }

    pub fn remove_genomes(& mut self, filename: &str){
        let file = File::open(filename).expect("ERROR: CAN NOT READ FILE\n");
        let reader = BufReader::new(file);
        let mut file_genomes: HashSet<String> = HashSet::new();
        for line in reader.lines() {
            let l = line.unwrap();
            file_genomes.insert(l);
        }
        let names = file_genomes.clone();


        let names_r:HashSet<String> = self.column_name.values().cloned().collect();
        let kk: HashSet<_> = names_r.difference(&names).collect();


        let mut to_remove = Vec::new();
        let mut rr = Vec::new();
        for x in self.column_name.iter(){
            if names.contains(x.1){
                to_remove.push(x.0.clone());
            }
        }

        for x in self.column_name.iter(){
            if kk.contains(x.1){
                rr.push(x.0.clone());
            }
        }

        println!("{:?}", to_remove);
        println!("this is k {:?}", rr);
        for x in rr.iter(){
            self.column_name.remove(&x);
        }

        rr.sort();
        println!("{:?}", rr);
        for (index, x) in rr.iter().enumerate(){
            self.matrix.matrix_core.remove((x.clone() as usize) - index);
        }
        println!("{:?}", rr);
    }

    pub fn make_binary(& mut self, thresh: u32){
        self.matrix_bin = self.matrix.copy(thresh);
    }

    pub fn write_bed(&self, out_prefix: &str, t: &str){
        //hexdump -C test.bin
        // xxd -b file
        // xxd file
        // SNP: 00000001 , 0
        // IND: 00000000, 1

        let mut buff: Vec<u8> = vec![108, 27, 1];
        let h2 = trans2( &self.matrix_bin);
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

    pub fn filter(&self) -> Vec<u32>{
        eprintln!("Filtering");
        println!("{} {}", self.matrix_bin.len(), self.matrix_bin[0].len());
        let k:Vec<Vec<bool>>= trans2(&self.matrix_bin);
        let mut k2 = Vec::new();
        let mut count = 0;

        let mut to_remove: Vec<u32> = Vec::new();

        for (i, x) in k.iter().enumerate(){
            let mut sum = 0;
            for y in x.iter() {
                if *y == true {
                    sum += 1;
                }
            }
            if (sum as usize != x.len()) | (sum == 0) {
                k2.push(x.clone());
            } else {
                println!("{} {}", sum, x.len());
                to_remove.push(i as u32);
                count += 1;
            }


        }
        eprintln!("Boring SNPs {}", to_remove.len());
        let k3 = trans2(&k2);

        eprintln!("Before {}  After {}", k[0].len(), k3[0].len());
        return to_remove;
    }


    /// Reduce binary shit
    pub fn reduce_comb1(& mut self) -> (Vec<usize>, Vec<usize>){
        // Meta
        // h1,h2 -> meta
        // hm -> BiMap (vec -> usize)
        let mut hm: BiMap<_,_> = BiMap::new();
        let mut h1: Vec<usize> = Vec::new();
        let mut h2: Vec<usize> = Vec::new();


        // Make SNPs Vector
        println!("len1 {:?}", self.matrix_bin.len());
        let k: Vec<Vec<bool>>= trans2(&self.matrix_bin);

        let mut count = 0;
        for (index, x) in k.iter().enumerate() {
            if ! hm.contains_left(x) {
                hm.insert(x.clone(), count);
                h1.push(index);
                h2.push(count);
                count += 1;
            } else {

                h1.push(index);
                h2.push(hm.get_by_left(x).unwrap().clone());
            }
        }
        let mut h : Vec<Vec<bool>> = Vec::new();
        for x in 0..hm.iter().len(){
            h.push(hm.get_by_right(&x).unwrap().clone());
        }


        println!("Reduce10 {}", count);
        println!("Reduce10 {:?}", hm);
        println!("Reduce101231 {:?}", h);

        self.matrix_bin = h;
        (h1,h2)
    }


    /// Write the names - helper function
    pub fn write_names(&self) {
        let f = File::create("holyshit1").expect("Unable to create file");
        let mut f = BufWriter::new(f);
        for x in self.column_name.iter(){
            write!(f, "{}\n", x.1).expect("Can not write file");
        }
    }



}


pub fn matrix_node10(gwrapper: &GraphWrapper, graph: &NGfa, mw: & mut MatrixWrapper2, h2: & mut BiMap<u32, usize>) {
    let mut h: Vec<u32> = graph.nodes.keys().cloned().collect();

    h.sort_by(|a, b| a.partial_cmp(b).unwrap());

    for (i, x) in h.iter().enumerate() {
        let node: u32 = x.clone();
        h2.insert(node, i);
    }

    for (index, (name, paths)) in gwrapper.genomes.iter().enumerate(){
        eprintln!("{}", name);
        mw.column_name.insert( index as u32, name.clone());
        let mut nody: Vec<u32> = vec![0; h2.len()] ;
        for x in paths.iter(){
            for y in x.nodes.iter(){
                nody[*h2.get_by_left(&y).unwrap()] += 1;

            }
        }
        mw.matrix.matrix_core.push(nody);
    }
}



pub fn matrix_node11(gwrapper: &GraphWrapper, graph: &NGfa, mw: & mut MatrixWrapper2, h2: & mut BiMap<u16, usize>) {
    let mut h: Vec<u32> = graph.nodes.keys().cloned().collect();

    h.sort_by(|a, b| a.partial_cmp(b).unwrap());

    for (i, x) in h.iter().enumerate() {
        let node: u16 = (*x as u16).clone();
        h2.insert(node, i);
    }

    for (index, (name, paths)) in gwrapper.genomes.iter().enumerate(){
        eprintln!("{}", name);
        mw.column_name.insert( index as u32, name.clone());
        let mut nody: Vec<u32> = vec![0; h2.len()] ;
        for x in paths.iter(){
            for y in x.nodes.iter(){
                nody[*h2.get_by_left(&(*y as u16) ).unwrap()] += 1;

            }
        }
        mw.matrix.matrix_core.push(nody);
    }

}

/// Make matrix for directed nodes // check this
pub fn matrix_dir_node2(gwrapper: &GraphWrapper, graph: &NGfa, mw: & mut MatrixWrapper2, j: & mut BiMap<(u32, bool), usize>){


    // Check all dir nodes
    let mut v: HashSet<(u32, bool)> = HashSet::new();
    for x in graph.paths.iter(){
        for x2 in 0..x.nodes.len(){
            v.insert((x.nodes[x2], x.dir[x2]));
        }
    }

    // Sort dirNod by node 1
    let mut k: Vec<_> = v.into_iter().collect();
    k.sort_by_key(|k| k.0);
    let mut count = 0;
    for x in k.iter(){
        j.insert(x.clone(), count);
        count += 1;
    }


    //eprintln!("{:?}", bim);







    for (index, (name, vec1)) in gwrapper.genomes.iter().enumerate(){
        mw.column_name.insert(index as u32, name.clone());
        //println!("{}", name);
        let mut dir_nodes : Vec<u32> = vec![0; j.len()] ;
        for x in vec1.iter(){
            for x2 in 0..x.nodes.len(){
                let node: u32 = x.nodes[x2].clone();

                dir_nodes[*j.get_by_left(&(node, x.dir[x2].clone())).unwrap()] += 1;

            }
        }
        mw.matrix.matrix_core.push(dir_nodes);
    }
}


/// Make matrix for edges
pub fn matrix_edge2(gwrapper: &GraphWrapper, graph: &NGfa, mw: & mut MatrixWrapper2, h2: & mut BiMap<(u32, bool, u32, bool), usize>){



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
        h2.insert(x.clone(), count);
        count += 1;
    }






    for (index, (name, vec1)) in gwrapper.genomes.iter().enumerate(){
        mw.column_name.insert( index as u32, name.clone());
        let mut dir_nodes : Vec<u32> = vec![0; h2.len()] ;
        for x in vec1.iter(){
            for x2 in 0..x.nodes.len()-1{
                let node: u32 = x.nodes[x2];
                let node2: u32 = x.nodes[x2+1];
                dir_nodes[*h2.get_by_left(&(node, x.dir[x2], node2, x.dir[x2+1])).unwrap()] += 1;

            }
        }
        mw.matrix.matrix_core.push(dir_nodes);
    }
}



/// Writing bim file
/// Information: https://www.cog-genomics.org/plink/1.9/formats#bim
pub fn write_bim<T>(ll: &BiMap<T, usize>, out_prefix: &str, t: &str)
where
T: Debug + std::hash::Hash + std::cmp::Eq + Ord
{
    let f = File::create([out_prefix, t,  "bim"].join(".")).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    let mut helper_vec = Vec::new();

    for (k,v) in ll.iter(){
        helper_vec.push((v,k));
    }
    helper_vec.sort_by_key(|k| k.1);



    for (k,_v) in helper_vec.iter(){
        write!(f, "{}\t{}\t{}\t{}\t{}\t{}\n", "graph", ".", 0, k, "A", "T").expect("Not able to write ");
    }






    let f = File::create([out_prefix, t,  "bimhelper"].join(".")).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    for (v,k) in helper_vec.iter(){
        write!(f, "{}\t{:?}\n", v, k).expect("Not able to write");

    }
}



/// Writing file wrapper
pub fn write_matrix(se: & mut MatrixWrapper2, what: &str, out_prefix: &str, t: &str){
    if (what == "bed") | (what == "all"){
        if se.matrix_bin.is_empty(){
            se.matrix_bin = se.matrix.copy(1);
        };
        se.write_bed(out_prefix, t);
        //se.write_bim(out_prefix, t);
        write_helper2(se, out_prefix);
    }
    if (what == "bimbam") | (what == "all"){
        //se.matrix.write_bimbam(out_prefix, t);
    }
}


pub fn write_helper2(se: & mut MatrixWrapper2, out_prefix: &str){
    let f = File::create([out_prefix,  "bimbim"].join(".")).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    for x in 0..se.column_name.len(){
        write!(f, "{}\n", se.column_name.get(&(x as u32)).unwrap()).expect("Can not write file");
    }
}



/// new matrix
/// Make matrix from bit vector
pub fn matrix_pack_bit(filename: &str, mw: & mut MatrixWrapper2, h2: & mut BiMap<u32, usize>) {
    let g: Vec<u8> = get_file_as_byte_vec(filename);
    let k: Vec<ReaderBit> = wrapper_bool(&g);

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
        h2.insert(x as u32, x);
    }
}



pub fn matrix_pack_u16(filename: &str, mw: & mut MatrixWrapper2, h2: & mut BiMap<u32, usize>) {
    let g: Vec<u8> = get_file_as_byte_vec(filename);
    let k: Vec<ReaderU16> = wrapper_u16(&g);

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
        h2.insert(x as u32, x);
    }
}

pub fn remove_bimap<T>(bm: & mut BiMap<T, usize>, v: Vec<u32>)
where
T:  Debug + std::hash::Hash + std::cmp::Eq
{

    for x in v.iter(){
        bm.remove_by_right(&(*x as usize));
    }

}
