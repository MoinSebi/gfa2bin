use bimap::BiMap;
use gfaR_wrapper::{GraphWrapper, NGfa};
use hashbrown::HashSet;
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
pub fn matrix_node(gwrapper: &GraphWrapper, graph: &NGfa, mw: & mut MatrixWrapper, h2: & mut BiMap<u32, usize>) {
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