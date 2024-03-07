
use crate::core::helper::{merge_u32_to_u64, to_string1, Feature};
use crate::r#mod::input_data::{FileData};
use clap::ArgMatches;
use gfa_reader::{NCGfa};
use hashbrown::HashSet;

use std::cmp::max;
use std::fs::File;
use std::io::Write;

pub fn find_main(matches: &ArgMatches) {
    let graph_file = matches.value_of("gfa").unwrap();
    let feature_file = matches.value_of("features").unwrap();
    let output = matches.value_of("output").unwrap();
    let length = matches.value_of("length").unwrap().parse::<i128>().unwrap();
    let data = FileData::from_file(feature_file);
    let feature = data.feature;

    // Read the graph
    let mut graph: NCGfa<()> = NCGfa::new();
    graph.parse_gfa_file(graph_file, false);
    let paths = &graph.paths;

    // Get the node size
    let node_size = node_size(&graph);

    //
    let mut position_nodesize = Vec::new();
    let mut vec_res_u64 = Vec::new();

    for path in paths.iter() {
        let mut vec_u64 = Vec::new();
        let mut index = Vec::new();
        let mut pos = 0;
        for i in 0..path.nodes.len() - 1 {
            index.push([pos, node_size[path.nodes[i] as usize]]);
            pos += node_size[path.nodes[i] as usize];
            let v1 = path.nodes[i];
            let v2 = path.dir[i];
            let v3 = path.nodes[i + 1];
            let v4 = path.dir[i + 1];

            if feature == Feature::Node {
                vec_u64.push(v1 as u64);
            } else if feature == Feature::DirNode {
                vec_u64.push(v1 as u64 * 2 + v2 as u64);
            } else if feature == Feature::Edge {
                let u1 = v1 * 2 + v2 as u32;
                let u2 = v3 * 2 + v4 as u32;
                vec_u64.push(merge_u32_to_u64(u1, u2));
            }
        }
        vec_res_u64.push(vec_u64);
        position_nodesize.push(index)
    }

    let file = File::create(output).unwrap();
    let mut writer = std::io::BufWriter::new(file);
    let data_hs = data.data.iter().collect::<HashSet<&u64>>();
    for (i, x) in vec_res_u64.iter().enumerate() {
        for (i2, y) in x.iter().enumerate() {
            if data_hs.contains(y) {
                writeln!(
                    writer,
                    "{}\t{}\t{}\tID:{};NS:{};NB:{}",
                    graph.paths[i].name,
                    max(0, position_nodesize[i][i2][0] as i128 - length),
                    position_nodesize[i][i2][0] as i128 + length,
                    to_string1(*y, &feature),
                    position_nodesize[i][i2][1],
                    position_nodesize[i][i2][0],
                )
                .expect("Error writing to file")
            }
        }
    }
}

// Size of each node
pub fn node_size(graph: &NCGfa<()>) -> Vec<usize> {
    let res = graph
        .nodes
        .iter()
        .map(|x| x.seq.len())
        .collect::<Vec<usize>>();
    res
}
