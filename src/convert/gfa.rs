use std::sync::{Arc, Mutex, MutexGuard};
use std::thread;
use bifurcation::helper::{chunk_inplace};
use bimap::BiMap;
use gfaR_wrapper::{GraphWrapper, NGfa, NPath};
use hashbrown::{HashMap, HashSet};
use log::{debug};
use crate::MatrixWrapper;


/// Multithread wrapper to get MatrixWrapper from a graph
pub fn matrix_node_wrapper(_gwrapper: & GraphWrapper, graph: &NGfa, mw: & mut MatrixWrapper, bimap: & mut Vec<u32>, threads: &usize) {
    // Create a VECTOR of all node id
    let mut nodes: Vec<u32> = graph.nodes.keys().cloned().collect();

    // Sort the nodes
    nodes.sort_by(|a, b| a.partial_cmp(b).unwrap());


    // Create BiMap of nodeID <> index
    let mut gz = HashMap::new();
    for (index, id) in nodes.iter().enumerate() {
        gz.insert(id.clone(), index);
    }

    // Chunking graph data for multithreading
    let chunks = chunk_inplace(graph.paths.clone(), threads.clone());
    // Add handles, result data structure, and the BiMap
    let mut handles = Vec::new();
    let result2 = Arc::new(Mutex::new(vec![]));
    let k2 = Arc::new(gz.clone());

    // Iterate over each chunk
    // Copy the result
    // Copy the BiMap
    for chunk in chunks{
        let r2 = result2.clone();
        let rro = k2.clone();
        let handle = thread::spawn(move || {

            for pair in chunk.iter(){
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
    combine_chromosomes(o, _gwrapper, mw);




}


/// Helper function for matrix_node_wrapper
///
/// Counting number of nodes in a path
pub fn matrix_node(path: &NPath, h2: &Arc<HashMap<u32, usize>>) -> Vec<u32>{
    let mut count: Vec<u32> = vec![0; h2.len()];

    for x in path.nodes.iter(){
        count[*h2.get(&x).unwrap()] += 1;
    }
    count
}


/// Matrix constructor from nodes
pub fn matrix_dir_node(_gwrapper: & GraphWrapper, graph: &NGfa, mw: & mut MatrixWrapper, bimap: &mut Vec<(u32, bool)>, threads: &usize){

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
    k.sort_by_key(|e| e.0);

    let mut map1 = HashMap::new();
    for (index, x) in k.iter().enumerate(){
        map1.insert( x.clone(), index);
    }


    let mut handles = Vec::new();
    let result2 = Arc::new(Mutex::new(vec![]));
    let k2 = Arc::new(map1.clone());

    for chunk in chunks {
        let r2 = result2.clone();
        let rro = k2.clone();
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
    bimap = k;

    let o = result2.lock().unwrap();
    combine_chromosomes(o, _gwrapper, mw);

}

pub fn tt (path: &NPath, h2: &Arc<HashMap<(u32, bool), usize>>) -> Vec<u32> {
    let mut dir_nodes : Vec<u32> = vec![0; h2.len()] ;
    for x2 in 0..path.nodes.len(){
        let node: u32 = path.nodes[x2].clone();

        dir_nodes[*h2.get(&(node, path.dir[x2].clone())).unwrap()] += 1;

    }
    println!("dir nodes {:?}", dir_nodes);
    return dir_nodes



}

/// Matrix constructor from edges
pub fn matrix_edge(_gwrapper: & GraphWrapper, graph: &NGfa, mw: & mut MatrixWrapper, bimap: &mut Vec<(u32, bool, u32, bool)>, threads: &usize){



    let mut v: HashSet<(u32, bool, u32, bool)> = HashSet::new();
    for x in graph.edges.iter(){
        let j: u32 = x.to;
        let j2: u32 = x.from;
        v.insert((j2, x.from_dir, j, x.to_dir));
    }

    let mut k: Vec<_> = v.into_iter().collect();
    k.sort_by_key(|k| k.0);
    let mut bw = HashMap::new();
    for (index, x) in k.iter().enumerate(){
        bw.insert(x.clone(), index);
    }



    let mut handles = Vec::new();
    let result2 = Arc::new(Mutex::new(vec![]));
    let k2 = Arc::new(bw.clone());
    let chunks = chunk_inplace(graph.paths.clone(), threads.clone());


    for chunk in chunks {
        let r2 = result2.clone();
        let rro = k2.clone();
        let handle = thread::spawn(move || {
            for x in chunk.iter(){
                let g = matrix_edge_helper(x, & rro);
                let mut rr = r2.lock().unwrap();
                rr.push((x.name.clone(), g));
            }
        });
        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap()

    }

    bimap = k;
    let o = result2.lock().unwrap();
    combine_chromosomes(o, _gwrapper, mw);


}

/// Helper function for matrix_edge
///
/// Counting number of edges in a path
pub fn matrix_edge_helper(path: &NPath, h2: &Arc<HashMap<(u32, bool, u32, bool), usize>>) -> Vec<u32> {
    let mut dir_nodes : Vec<u32> = vec![0; h2.len()] ;

    for x2 in 0..path.nodes.len()-1{
        let node: u32 = path.nodes[x2];
        let node2: u32 = path.nodes[x2+1];
        dir_nodes[*h2.get(&(node, path.dir[x2], node2, path.dir[x2+1])).unwrap()] += 1;

    }
    dir_nodes
}

/// Combine path of the same genome
///
/// Based on genome name in PanSN naming
pub fn combine_chromosomes(o: MutexGuard<Vec<(String, Vec<u32>)>>, h: & GraphWrapper, mw: &mut MatrixWrapper){
    let mut data1 = HashMap::new();
    for (pathname, counts) in o.iter(){
        let genome_name = h.path2genome.get(pathname).unwrap();
        data1.entry(genome_name).or_insert(vec![counts.clone()]).push(counts.clone());
    }

    let mut data_merge = HashMap::new();
    for (key, value) in data1.iter_mut(){
        let mut f = value[0].clone();
        for x in value.iter().skip(1){
            f = f.iter().zip(x.iter()).map(|(x, y)| x + y).collect();
        }
        data_merge.insert(key, f);
    }

    mw.matrix_core = Vec::with_capacity(o.len());
    for (i, x) in data_merge.iter().enumerate(){
        mw.matrix_core.push(x.1.clone());
        mw.column_name.insert(i as u32, (**x.0).clone());
    }

}

