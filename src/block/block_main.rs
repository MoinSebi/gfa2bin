use crate::core::core::MatrixWrapper;
use crate::core::helper::{merge_u32_to_u64, Feature};
use crate::subpath::subpath_main::traversal2samples;
use clap::ArgMatches;
use gfa_reader::{Gfa, Pansn};

use log::info;
use std::collections::HashSet;
use std::fmt::Error;
use std::{fmt, io};

/// Block main function
///
/// Easy block function
/// Extract the subpath from a graph for each node
pub fn block_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // Input
    let graph_file = matches.value_of("gfa").unwrap();
    let sep = matches.value_of("PanSN").unwrap_or(" ");

    // Block parameters
    let window: usize = matches
        .value_of("window")
        .unwrap()
        .parse::<usize>()
        .unwrap();
    let step_size: usize = matches.value_of("step").unwrap().parse().unwrap();
    let sequence = matches.value_of("sequence");
    let sequence_window = matches.value_of("sequence window");

    let cutoff_distance: usize = matches.value_of("distance").unwrap().parse().unwrap();

    // Output
    let output_prefix = matches.value_of("output").unwrap();
    let split = matches
        .value_of("split")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();

    info!("Block subcommand");
    info!("Graph file: {}", graph_file);
    info!("Separator: {}", sep);
    info!("Window size: {}", window);
    info!("Distance: {}", cutoff_distance);
    info!("Step size: {}", step_size);
    info!("Split size: {}", split);
    info!("Output prefix: {}\n", output_prefix);

    info!("Reading graph file");
    let window = window;
    let graph: Gfa<u32, (), ()> = Gfa::parse_gfa_file(graph_file);
    let wrapper: Pansn<u32, (), ()> = Pansn::from_graph(&graph.paths, sep);

    info!("Indexing graph");
    let a = block_wrapper(&graph, step_size, window, sequence, sequence_window)?;
    info!("Number of blocks: {}", a.len());

    let b = node_size(&graph);

    info!("Extracting blocks");
    let mw = wrapper_blocks(&wrapper, b, a, cutoff_distance);

    info!("Writing blocks");
    let chunk_size = (mw.matrix_bit.len() / split) + 1;
    let chunks = mw.matrix_bit.chunks(chunk_size);

    let len = chunks.len();
    for (index, _y) in chunks.enumerate() {
        mw.write_fam(index, output_prefix, mw.feature, len, 0.0);
        mw.write_bed(index, output_prefix, mw.feature, len);
        mw.write_bim(index, output_prefix, &mw.feature, len);
    }
    Ok(())
}


pub fn block_wrapper(graph: &Gfa<u32, (), ()>, step: usize, window_size: usize, sequence: Option<&str>, sequence_window: Option<&str>) -> Result<Vec<[u32; 2]>, Box<dyn std::error::Error>> {
    let mut blocks = Vec::new();
    if sequence.is_none() {
        blocks = blocks_node(graph, step, window_size);
    } else {
        let sequence = sequence.unwrap().parse()?;
        let sequence_window = sequence_window.unwrap().parse()?;
        blocks = block_seq(graph, sequence, sequence_window);
    }
    Ok(blocks)
}

/// Make blocks
///
///  - A block starts at a node and end at a node
///  - Returns start and end nodes of a block
pub fn blocks_node(graph: &Gfa<u32, (), ()>, step: usize, wsize: usize) -> Vec<[u32; 2]> {
    let mut gg = Vec::new();
    for x in (0..(graph.segments.len()-wsize)).step_by(step) {
        gg.push([graph.segments[x].id, graph.segments[x+wsize].id]);
    }
    gg
}


pub fn block_seq(graph: &Gfa<u32, (), ()>, size: usize, step_size: usize) -> Vec<[u32; 2]>{
    let starts = get_starts(graph, step_size);
    let mut result = Vec::new();
    let mut sequence = 0;
    let mut starting_id = graph.segments[0].id;
    for x in starts.iter() {
        for y in *x..graph.segments.len() {
            sequence += graph.segments[y].sequence.get_len();
            if sequence > size {
                result.push([starting_id, graph.segments[y].id]);
                starting_id = graph.segments[y].id;
                sequence = 0;
            }
        }
    }
    result
}

pub fn get_starts(graph: &Gfa<u32, (), ()>, step: usize) -> Vec<usize> {
    let mut gg = Vec::new();
    let mut pos = 0;
    for (index, x) in graph.segments.iter().enumerate() {
        pos += x.sequence.get_len();
        if pos > step {
            gg.push(index);
            pos = 0;
        }
    }
    gg
}



/// Node size index
pub fn node_size(graph: &Gfa<u32, (), ()>) -> Vec<usize> {
    let mut node_size = Vec::new();
    for node in graph.segments.iter() {
        node_size.push(node.sequence.get_len());
    }
    node_size
}

/// Wrapper function for blocks
///
/// Iterate over blocks
/// Iterate over each path
/// Iterate over each node
/// If node is in block:
///     if block is not empty: add to block
///    else: create new block
/// else: add distance
pub fn wrapper_blocks(
    graph2: &Pansn<u32, (), ()>,
    node_size: Vec<usize>,
    block: Vec<[u32; 2]>,
    max_distance: usize,
) -> MatrixWrapper {
    // this is the output
    let mut mw = MatrixWrapper::new();

    // Iterate over the blocks
    for x in block.iter() {
        let block_hashset = (x[0]..x[1]).collect::<HashSet<u32>>();
        let mut all_blocks = Vec::new();

        // Iterate over each
        for (genome_id, path) in graph2.genomes.iter().enumerate() {
            for (haplo_id, x1) in path.haplotypes.iter().enumerate() {
                for (path_id, x) in x1.paths.iter().enumerate() {
                    //
                    let mut block_array: [usize; 3] = [0; 3]; // Triple 0
                    let mut distance = 0;
                    for (i, node) in x.nodes.iter().enumerate() {
                        if block_hashset.contains(node) {
                            if block_array[2] == 0 {
                                block_array[0] = i;
                                block_array[1] = i;
                                block_array[2] = node_size[*node as usize - 1];
                                distance = 0;
                            } else {
                                block_array[1] = i;
                                block_array[2] += node_size[*node as usize - 1];
                                distance = 0;
                            }
                        } else {
                            distance += node_size[*node as usize - 1];
                            if distance > max_distance && block_array[2] != 0 {
                                all_blocks.push((
                                    genome_id,
                                    haplo_id,
                                    path.haplotypes.len(),
                                    &x.nodes[block_array[0]..block_array[1]],
                                ));
                                block_array = [0; 3];
                            }
                        }
                    }
                    if block_array[2] != 0 {
                        all_blocks.push((
                            genome_id,
                            haplo_id,
                            path.haplotypes.len(),
                            &x.nodes[block_array[0]..block_array[1]],
                        ));
                    }
                }
            }
        }

        //
        all_blocks.sort();
        if all_blocks.is_empty() {
            continue;
        }

        let bb = all_blocks.len();

        //mw.matrix_bit.extend(traversal2samples(all_blocks, 20));
        mw.geno_names
            .extend((0..bb).map(|aa| merge_u32_to_u64(x[0] / 2, aa as u32)));
    }
    mw.feature = Feature::Block;
    mw.window_size = (block[0][1] - block[0][0]) as usize;
    mw.shape = (mw.matrix_bit.len(), mw.matrix_bit[0].len());
    mw.sample_names = graph2.genomes.iter().map(|x| x.name.clone()).collect();
    mw
}
