use std::collections::HashMap;
use gfaR::{Gfa, Path};
use std::borrow::{Borrow, BorrowMut};
use std::fmt::{Debug, Display};
use std::collections::hash_set::Iter;
use std::fs::File;
use std::io::{Write, BufWriter};
use gfaR_wrapper::{NGfa, GraphWrapper};

#[derive(Debug, Clone)]
pub struct matrix{
    pub shape: (u32, u32),
    pub matrix_core: Vec<Vec<u32>>,
}

impl matrix {
    pub fn new() -> Self {
        let shape: (u32, u32) = (0,0);
        let matrix: Vec<Vec<u32>> = Vec::new();
        Self {
            shape: shape,
            matrix_core: matrix,
        }
    }

    pub fn normalize(& self) -> Vec<Vec<bool>>{
        let mut new_mat: Vec<Vec<bool>> = Vec::new();
        for x in self.matrix_core.iter(){
            let mut k: Vec<bool> = Vec::new();
            for y in x.iter(){
               k.push(y.clone() != 0);
            }
            new_mat.push(k);
        }
        new_mat
    }

    pub fn copy(&self, number: u32) -> Vec<Vec<bool>>{
        let mut new_mat: Vec<Vec<bool>> = Vec::new();
        for x in self.matrix_core.iter(){
            let mut k: Vec<bool> = Vec::new();
            for y in x.iter(){
                k.push(y.clone() >= 0);

            }
            new_mat.push(k);
        }
        new_mat
    }



    pub fn min_max(&self) -> Vec<Vec<f32>>{
        let mut new_mat:  Vec<Vec<f32>> = Vec::new();
        for x in 0..self.matrix_core[0].len(){
            let mut max = 0;
            for y in 0..self.matrix_core.len(){
                if self.matrix_core[y][x] > max{
                    max = self.matrix_core[y][x]
                }
            }
            let mut vec_new: Vec<f32> = Vec::new();
            for y in 0..self.matrix_core.len(){
                vec_new.push(self.matrix_core[y][x] as f32/max as f32) //this is wrong
            }
            new_mat.push(vec_new);
        }

        new_mat
    }


    pub fn write_bimbam(&self){
        println!("daskdhaskjhdjashdjas");
        let f = File::create("test.bimbam").expect("Unable to create file");
        let mut f = BufWriter::new(f);
        let k = self.min_max();
        for (index, x) in k.iter().enumerate(){
            let j: Vec<String> = x.iter().map(|i| format!("{}", i)).collect();
            write!(f, "{}\t{}\t{}\t", index, "A", "T");
            write!(f, "{}\n", j.join("\t"));
        }
    }

    /// Add buffer here
    pub fn write_bed(&self, out_prefix: &str){
        //hexdump -C test.bin
        // xxd -b file
        // xxd file
        // SNP: 00000001 , 0
        // IND: 00000000, 1
        let dummy_vec1: Vec<u8> = vec![108, 27, 1];

        let h = self.normalize();
        println!("{}", self.matrix_core.len());
        let j: Vec<&[bool]> = h[0].chunks(4).collect();
        let mut j2: Vec<u8> = vec![108, 27, 1];
        for x in j{
            j2.push(binary2dec_bed(x));
        }
        let mut file = File::create([out_prefix, "bed"].join(".")).expect("Hilfe hilfe");
        file.write_all(&j2).expect("hilfe2");

    }



}


pub fn binary2dec_bed(vecc: &[bool]) -> u8{
    let mut result: u8 = 0;
    let mut count = 0;
    for x in vecc.iter().rev(){
        let t: u8 = 2;
        result += (t.pow(count as u32)) * (*x as u8);
        count += 1;
        result += (t.pow(count as u32)) * (*x as u8);
        count += 1;
    }
    result
}







pub struct matrix_wrapper<T: Debug>{
    pub matrix: matrix,
    pub column_name: HashMap<String, u32>,
    pub row_name: HashMap<T, usize>,

}

impl <T>matrix_wrapper<T>
where
T: Debug
{
    pub fn new() -> Self{
        let matrx = matrix::new();
        let col: HashMap<String, u32> = HashMap::new();
        let row: HashMap<T, usize> = HashMap::new();
        Self{
            matrix: matrx,
            column_name: col,
            row_name: row,
        }
    }

    pub fn bim(&self, out_prefix: &str){
        let f = File::create("test.bim").expect("Unable to create file");
        let mut f = BufWriter::new(f);
        for (k,v) in self.row_name.iter(){
            //write!(f, "{}\t{}\t{:?}\t{}\t{}\n", "graph", 0, k, v, v)
            write!(f, "{}\t{}\t{:?}\t{}\t{}\n", 0, 0, 0,0,0);
        }

    }

    pub fn bim_helper(&self){
        let f = File::create("test.bimhelper").expect("Unable to create file");
        let mut f = BufWriter::new(f);
        for (k,v) in self.row_name.iter(){
            write!(f, "{}\t{}\t{:?}\t{}\t{}\n", "graph", 0, k, v, v);
        }
    }




}


// pub fn test11(graph: &NGfa) -> matrix_wrapper<u32>{
//     let mut mw: matrix_wrapper<u32> = matrix_wrapper::new();
//
//
//     for (i, x) in graph.nodes.iter().enumerate(){
//         mw.row_name.insert (x.0.clone(), i);
//     }
//     println!("Len is {}", mw.row_name.len());
//     for (index, path) in graph.paths.iter().enumerate(){
//         mw.column_name.insert(path.name.clone(), index as u32);
//         let mut nody: Vec<u32> = vec![0; mw.row_name.len()] ;
//         println!("{}", nody.len());
//         for x in path.nodes.iter(){
//             nody[mw.row_name[x]] += 1;
//
//         }
//         mw.matrix.matrix_core.push(nody);
//     }
//     mw
//
// }


pub fn test1(gwrapper: &GraphWrapper, graph: &NGfa) -> matrix_wrapper<u32>{
    let mut mw: matrix_wrapper<u32> = matrix_wrapper::new();


    for (i, x) in graph.nodes.iter().enumerate(){
        let node:u32 = x.0.clone();
        mw.row_name.insert (node, i);
    }


    for (index, (name, paths)) in gwrapper.genomes.iter().enumerate(){

        mw.column_name.insert(name.clone(), index as u32);
        let mut nody: Vec<u32> = vec![0; mw.row_name.len()] ;
        println!("{}", nody.len());
        for x in paths.iter(){
            for y in x.nodes.iter(){
                nody[mw.row_name[&y]] += 1;

            }
        }
        mw.matrix.matrix_core.push(nody);
    }
    mw

}

pub fn test22(gwrapper: &GraphWrapper, graph: &NGfa) -> matrix_wrapper<(u32, bool)>{
    let mut mw: matrix_wrapper<(u32, bool)> = matrix_wrapper::new();
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
    //println!("{:?}", mw.row_name);



    mw


}

pub fn test2(graph: &Gfa) -> matrix_wrapper<(u32, bool)>{
    let mut mw: matrix_wrapper<(u32, bool)> = matrix_wrapper::new();

    let mut count = 0;
    for x in graph.paths.iter(){
        for x2 in 0..x.nodes.len(){
            let node: u32 = x.nodes[x2].parse().unwrap();
            if ! mw.row_name.contains_key(&(node, x.dir[x2])){
                mw.row_name.insert((node, x.dir[x2].clone()), count);
            }
        }
    }

    for (index, x) in graph.paths.iter().enumerate(){
        mw.column_name.insert(x.name.clone(), index as u32);
        let mut dir_nodes : Vec<u32> = vec![0; mw.row_name.len()] ;
        for x2 in 0..x.nodes.len(){
            let node: u32 = x.nodes[x2].parse().unwrap();
            dir_nodes[mw.row_name[&(node, x.dir[x2].clone())]] += 1;
        }
    }

    mw
}

pub fn test3(graph: &Gfa) -> matrix_wrapper<(u32, bool, u32, bool)>{
    let mut mw: matrix_wrapper<(u32, bool, u32, bool)> = matrix_wrapper::new();

    for (i, x) in graph.edges.iter().enumerate() {
        let j: u32 = x.to.parse().unwrap();
        let j2: u32 = x.from.parse().unwrap();
        mw.row_name.insert((j2, x.from_dir, j, x.to_dir), i);
    }

    for (index, x)  in graph.paths.iter().enumerate(){
        mw.column_name.insert(x.name.clone(), index as u32);
        let mut dir_nodes : Vec<u32> = vec![0; mw.row_name.len()] ;
        for x2 in 0..x.nodes.len()-1{
            let node: u32 = x.nodes[x2].parse().unwrap();
            let node2: u32 = x.nodes[x2+1].parse().unwrap();
            dir_nodes[mw.row_name[&(node, x.dir[x2], node2, x.dir[x2+1])]] += 1;
        }
    }

    mw
}
