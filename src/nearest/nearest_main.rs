use std::fs::File;
use std::hash::Hash;
use clap::ArgMatches;
use gfa_reader::Gfa;
use hashbrown::{HashMap, HashSet};
use crate::core::core::MatrixWrapper;
use crate::core::helper::Feature;
use crate::remove::input_data::FileData;

pub fn nearest_main(matches: &ArgMatches) {
    let graph_file = matches.value_of("gfa").unwrap();
    let feature_file = matches.value_of("nodes").unwrap();
    let output = matches.value_of("output").unwrap();
    let reference = matches
        .value_of("reference")
        .unwrap();

}

pub fn get_ref_nodes(graph: &Gfa<u32, (), ()>, names: &Vec<String>) -> HashSet<u32>{
    let mut nodes_hashset = HashSet::new();
    for path in graph.paths.iter(){
        if names.contains(&path.name){
            nodes_hashset.extend(&path.nodes)
        }
    }
    nodes_hashset
}

pub fn init_hm(graph: &Gfa<u32, (), ()>, checked_nodes: &HashSet<u32>) -> HashMap<u32, [u32; 2]>{
    let mut result_hm = HashMap::with_capacity(graph.segments.len());
    if !checked_nodes.is_empty(){
        for x in graph.segments.iter(){
            result_hm.insert(x.id, [x.id, u32::MAX]);
        }
    } else {
        for x in checked_nodes.iter(){
            result_hm.insert(*x, [*x, u32::MAX]);
        }
    }
    result_hm
}

pub fn read_nodes(graph: &Gfa<u32, (), ()>, names: &Vec<String>, checked_nodes: &HashSet<u32>) -> HashMap<u32, [u32; 2]>{
    let mut distance = 0;
    let mut reference_node = 0;
    let nodes_hs = get_ref_nodes(graph, names);
    let mut result_hm = init_hm(graph);
    for path in graph.paths.iter(){
        if !names.contains(&path.name){
            for node in path.nodes{
                if nodes_hs.contains(&node){
                    distance = 0;
                    reference_node = node;
                    result_hm[node] = [node, 0]
                } else {
                    if checked_nodes.contains(&node) {
                        if result_hm[node][1] > distance {
                            result_hm[node] = [reference_node, distance]
                        }
                    }
                    distance += graph.get_node_by_id(&node).length;

                }
            }
        }
    }
    result_hm
}


pub fn write_file(result: HashMap<u32, [u32; 2]>, output: &str, checked_nodes: &Vec<u32>, graph: &Gfa<u32, (), ()>, names: &Vec<String>) -> Result<(), std::io::Error>{

    let file_out = File::create(output)?;
    let mut output_reader = std::io::BufWriter::new(file_out);
    let setter = setter1(graph, names);
    writeln!(
            output_reader,
            "{}\t{}\t{}\t{}\t{}",
            "node", "ref_node", "ref_fasta", "ref_pos", "distance"
        )?;
    for x in checked_nodes.iter(){
        writeln!(
            output_reader,
            "{}\t{}\t{}\t{}\t{}",
            x, result[x][0], get_string(&setter,  result[x][0]), result[x][1], result[x][1]
        )?;
    }

    Ok(())

}

pub fn setter1(graph: &Gfa<u32, (), ()>, names: &Vec<String>) ->  Vec<(String, HashSet<u32>)>{
    let mut vv: Vec<(String, HashSet<u32>)> = Vec::new();
    for x in graph.paths.iter(){
        if names.contains(&x.name){
            vv.push((x.name.clone(), x.nodes.clone().iter().collect::<HashSet<u32>>()));
        }
    }
    return vv
}

pub fn get_string(aa:  &Vec<(String, HashSet<u32>)>, node: u32) -> String{
    let mut res: Vec<&String> = Vec::new();
    for x in aa.iter(){
        if x.1.contains(&node){
            res.push(&x.0);
        }
    }
    res.join(",")
}

pub fn get_pos(graph: Gfa<u32, (), ()>){
    let mut a = HashMap::new();
    for path in graph.paths.iter(){
        let mut o: HashMap<&u32, Vec<u32>> = HashMap::new();
        for p in path.nodes.iter(){
            if o.contains_key(p){
                let o = o.entry(p);
            } else {
                o.insert(p, vec![0]);
            }
        }
        a.
    }

}