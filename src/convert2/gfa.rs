use std::sync::{Arc, Mutex};
use std::thread;
use bifurcation::helper::{chunk_inplace, get_all_pairs};
use bimap::BiMap;
use bitvec::ptr::Mut;
use gfaR_wrapper::{GraphWrapper, NGfa, NPath};
use hashbrown::HashSet;
use log::debug;
use crate::MatrixWrapper;


impl MatrixWrapper{
    pub fn matrix_node(gwrapper: &GraphWrapper, graph: &NGfa, mm: & mut Self, h2: & mut BiMap<u32, usize>) {
        let mut h: Vec<u32> = graph.nodes.keys().cloned().collect();

        h.sort_by(|a, b| a.partial_cmp(b).unwrap());

        for (i, x) in h.iter().enumerate() {
            let node: u32 = x.clone();
            h2.insert(node, i);
        }

        for (index, (name, paths)) in gwrapper.genomes.iter().enumerate(){
            eprintln!("{}", name);
            mm.column_name.insert( index as u32, name.clone());
            let mut nody: Vec<u32> = vec![0; h2.len()] ;
            for x in paths.iter(){
                for y in x.nodes.iter(){
                    nody[*h2.get_by_left(&y).unwrap()] += 1;

                }
            }
            mm.matrix.matrix_core.push(nody);
        }
    }
}

/// Matrix constructor from nodes
/// # Arguments
/// * 'gwrapper' - Graph wrapper data structure
pub fn matrix_node_wrapper(gwrapper: &GraphWrapper, graph: &NGfa, mw: & mut MatrixWrapper, h2: & mut BiMap<u32, usize>) {
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


pub fn matrix_node_wrapper2(gwrapper: &GraphWrapper, graph: &NGfa, mw: & mut MatrixWrapper, h2: & mut BiMap<u32, usize>, threads: &usize) {
    let mut h: Vec<u32> = graph.nodes.keys().cloned().collect();

    h.sort_by(|a, b| a.partial_cmp(b).unwrap());

    for (i, x) in h.iter().enumerate() {
        let node: u32 = x.clone();
        h2.insert(node, i);
    }

    let chunks = chunk_inplace(graph.paths.clone(), threads.clone());


    let result = Arc::new(Mutex::new(mw.clone()));
    let wrapper = Arc::new(gwrapper);
    let mut handles = Vec::new();
    let result2 = Arc::new(Mutex::new(vec![]));
    let k = Arc::new(h2.clone());


    for chunk in chunks{
        let r = result.clone();
        let r2 = result2.clone();
        let wra = wrapper.clone();
        let rro = k.clone();
        let handle = thread::spawn(move || {
            for pair in chunk.iter(){
                debug!("Pair: {} ", pair.name);

                let mut h = matrix_node(pair, &rro);
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

    for (i, x) in o.iter().enumerate(){
        mw.matrix.matrix_core.push(x.1.clone());
        mw.column_name.insert(i as u32, x.0.clone());
    }




}

pub fn matrix_node(path: &NPath, h2: &Arc<BiMap<u32, usize>>) -> Vec<u32>{
    let mut nody: Vec<u32> = vec![0; h2.len()];
    for x in path.nodes.iter(){
        nody[*h2.get_by_left(&x).unwrap()] += 1;


    }
    nody
}



/// Matrix constructor from nodes
pub fn matrix_dir_node(gwrapper: &GraphWrapper, graph: &NGfa, mw: & mut MatrixWrapper, j: & mut BiMap<(u32, bool), usize>){



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


/// Matrix constructor from edges
pub fn matrix_edge(gwrapper: &GraphWrapper, graph: &NGfa, mw: & mut MatrixWrapper, h2: & mut BiMap<(u32, bool, u32, bool), usize>){



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