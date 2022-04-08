use std::sync::{Arc, Mutex};
use std::thread;
use bifurcation::helper::{chunk_inplace};
use bimap::BiMap;
use gfaR_wrapper::{NGfa, NPath};
use hashbrown::{HashSet};
use log::{debug};
use crate::MatrixWrapper;


/// Multithread wrapper to get MatrixWrapper from a graph
pub fn matrix_node_wrapper(graph: &NGfa, mw: & mut MatrixWrapper, bimap: & mut BiMap<u32, usize>, threads: &usize) {
    // Get all the nodes
    let mut h: Vec<u32> = graph.nodes.keys().cloned().collect();

    // Sort the nodes
    h.sort_by(|a, b| a.partial_cmp(b).unwrap());


    for (i, x) in h.iter().enumerate() {
        let node: u32 = x.clone();
        bimap.insert(node, i);
    }
    println!("path {:?}", graph.paths);

    // This is for multithreading
    let chunks = chunk_inplace(graph.paths.clone(), threads.clone());

    let mut handles = Vec::new();
    let result2 = Arc::new(Mutex::new(vec![]));
    let k = Arc::new(bimap.clone());
    println!("chunks {}", chunks.len());


    for chunk in chunks{
        let r2 = result2.clone();
        let rro = k.clone();
        let handle = thread::spawn(move || {

            for pair in chunk.iter(){
                debug!("Pair: {} ", pair.name);
                println!("djkalsjdla {}", pair.name);
                let h = matrix_node(pair, &rro);
                let mut rr = r2.lock().unwrap();
                rr.push((pair.name.clone(), h));

            }
        });
        handles.push(handle);

    }

    for handle in handles {
        handle.join().unwrap()

    }

    let o = result2.lock().unwrap();

    mw.matrix.matrix_core = Vec::with_capacity(o.len());
    println!("{:?}", mw.matrix.matrix_core);
    for (i, x) in o.iter().enumerate(){
        println!("dasdasd {:?}", x);
        mw.matrix.matrix_core.push(x.1.clone());
        mw.column_name.insert(i as u32, x.0.clone());
    }

    println!("{:?}", mw.matrix.matrix_core);




}

pub fn matrix_node(path: &NPath, h2: &Arc<BiMap<u32, usize>>) -> Vec<u32>{
    let mut count: Vec<u32> = vec![0; h2.len()];

    for x in path.nodes.iter(){
        count[*h2.get_by_left(&x).unwrap()] += 1;


    }
    count
}
/// Matrix constructor from nodes
pub fn matrix_dir_node(graph: &NGfa, mw: & mut MatrixWrapper, bimap: & mut BiMap<(u32, bool), usize>, threads: &usize){

    let chunks = chunk_inplace(graph.paths.clone(), threads.clone());

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
        bimap.insert(x.clone(), count);
        count += 1;
    }


    let mut handles = Vec::new();
    let result2 = Arc::new(Mutex::new(vec![]));
    let k = Arc::new(bimap.clone());

    for chunk in chunks {
        let r2 = result2.clone();
        let rro = k.clone();
        let handle = thread::spawn(move || {
            for x in chunk.iter(){
                let g = tt(x, & rro);
                let mut rr = r2.lock().unwrap();
                rr.push((x.name.clone(), g));
            }
        });
        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap()

    }

    let o = result2.lock().unwrap();

    mw.matrix.matrix_core = Vec::with_capacity(o.len());
    println!("len is {}", o.len());
    for (i, x) in o.iter().enumerate(){
        mw.matrix.matrix_core.push(x.1.clone());
        mw.column_name.insert(i as u32, x.0.clone());
    }
}

pub fn tt (path: &NPath, h2: &Arc<BiMap<(u32, bool), usize>>) -> Vec<u32> {
    let mut dir_nodes : Vec<u32> = vec![0; h2.len()] ;
    for x2 in 0..path.nodes.len(){
        let node: u32 = path.nodes[x2].clone();

        dir_nodes[*h2.get_by_left(&(node, path.dir[x2].clone())).unwrap()] += 1;

    }
    println!("{:?}", dir_nodes);
    return dir_nodes



}

/// Matrix constructor from edges
pub fn matrix_edge(graph: &NGfa, mw: & mut MatrixWrapper, bimap: & mut BiMap<(u32, bool, u32, bool), usize>, threads: &usize){



    let mut v: HashSet<(u32, bool, u32, bool)> = HashSet::new();
    for x in graph.edges.iter(){
        let j: u32 = x.to;
        let j2: u32 = x.from;
        v.insert((j2, x.from_dir, j, x.to_dir));
    }

    let mut k: Vec<_> = v.into_iter().collect();
    k.sort_by_key(|k| k.0);
    let mut count = 0;
    for x in k.iter(){
        bimap.insert(x.clone(), count);
        count += 1;
    }



    let mut handles = Vec::new();
    let result2 = Arc::new(Mutex::new(vec![]));
    let k = Arc::new(bimap.clone());
    let chunks = chunk_inplace(graph.paths.clone(), threads.clone());


    for chunk in chunks {
        let r2 = result2.clone();
        let rro = k.clone();
        let handle = thread::spawn(move || {
            for x in chunk.iter(){
                let g = tt2(x, & rro);
                let mut rr = r2.lock().unwrap();
                rr.push((x.name.clone(), g));
            }
        });
        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap()

    }

    let o = result2.lock().unwrap();

    mw.matrix.matrix_core = Vec::with_capacity(o.len());
    for (i, x) in o.iter().enumerate(){
        mw.matrix.matrix_core.push(x.1.clone());
        mw.column_name.insert(i as u32, x.0.clone());
    }

}

pub fn tt2(path: &NPath, h2: &Arc<BiMap<(u32, bool, u32, bool), usize>>) -> Vec<u32> {
    let mut dir_nodes : Vec<u32> = vec![0; h2.len()] ;

    for x2 in 0..path.nodes.len()-1{
        let node: u32 = path.nodes[x2];
        let node2: u32 = path.nodes[x2+1];
        dir_nodes[*h2.get_by_left(&(node, path.dir[x2], node2, path.dir[x2+1])).unwrap()] += 1;

    }
    dir_nodes
}

