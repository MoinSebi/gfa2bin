use crate::core::core::MatrixWrapper;
use crate::core::helper::{merge_u32_to_u64, Feature};
use crate::subpath::subpath_main::{function1, gfa_index, subpath_wrapper};
use clap::ArgMatches;
use gfa_reader::{NCGfa, NCPath, Pansn};
use hashbrown::HashMap;
use log::info;
use std::collections::HashSet;
use std::ffi::c_ushort;

/// Block main function
///
/// Easy block function
/// Extract the subpath from a graph for each node
pub fn block_main(matches: &ArgMatches) {

    // Input
    let graph_file = matches.value_of("gfa").unwrap();
    let sep = matches.value_of("PanSN").unwrap_or(" ");

    // Block param,eters
    let window: usize = matches
        .value_of("window")
        .unwrap()
        .parse::<usize>()
        .unwrap();
    let steps: usize = matches.value_of("step").unwrap().parse().unwrap();
    let distance: usize = matches.value_of("distance").unwrap().parse().unwrap();

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
    info!("Distance: {}", distance);
    info!("Step size: {}", steps);
    info!("Split size: {}", split);
    info!("Output prefix: {}\n", output_prefix);

    info!("Reading graph file");
    let window = window / 2;
    let mut graph: NCGfa<()> = NCGfa::new();
    graph.parse_gfa_file(graph_file, false);
    let wrapper: Pansn<NCPath> = Pansn::from_graph(&graph.paths, sep);

    info!("Indexing graph");
    let a = blocks_node(&graph, steps, window);
    info!("Number of blocks: {}", a.len());

    let b = node_size(&graph);

    info!("Extracting blocks");
    let mw = wrapper_blocks(&wrapper, b, a, distance);


    info!("Writing blocks");
    let chunk_size = (mw.matrix_bit.len() / split) + 1;
    let chunks = mw.matrix_bit.chunks(chunk_size);

    let len = chunks.len();
    for (index, _y) in chunks.enumerate() {
        mw.write_fam(index, output_prefix, mw.feature, len);
        mw.write_bed(index, output_prefix, mw.feature, len);
        mw.write_bim(index, output_prefix, &mw.feature, len);
    }
}

/// Make blocks
///
///  - A block starts at a node and end at a node
///  - Returns start and end nodes of a block
pub fn blocks_node(graph: &NCGfa<()>, step: usize, wsize: usize) -> Vec<[u32; 2]> {
    let glen = graph.nodes.len();
    let mut gg = Vec::new();
    for x in (wsize..glen - wsize).step_by(step) {
        gg.push([(x - wsize) as u32, (x + wsize) as u32]);
    }
    return gg;
}

/// Node size index
pub fn node_size(graph: &NCGfa<()>) -> Vec<usize> {
    let mut node_size = Vec::new();
    for node in graph.nodes.iter() {
        node_size.push(node.seq.len());
    }
    return node_size;
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
    graph2: &Pansn<NCPath>,
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
                    let mut block_array: [usize; 3] = [0; 3];       // Triple 0
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
                            if distance > max_distance {
                                if block_array[2] != 0 {
                                    all_blocks.push((
                                        genome_id,
                                        haplo_id,
                                        path_id,
                                        &x.nodes[block_array[0]..block_array[1]],
                                    ));
                                    block_array = [0; 3];
                                }
                            }
                        }
                    }
                    if block_array[2] != 0 {
                        all_blocks.push((
                            genome_id,
                            haplo_id,
                            path_id,
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

        mw.matrix_bit.extend(function1(all_blocks, 20));
        mw.geno_names
            .extend((0..bb).map(|aa| merge_u32_to_u64(x[0] / 2, aa as u32)));
    }
    mw.feature = Feature::Block;
    mw.window_size = (block[0][1] - block[0][0]) as usize;
    mw.shape = (mw.matrix_bit.len(), mw.matrix_bit[0].len());
    mw.sample_names = graph2.genomes.iter().map(|x| x.name.clone()).collect();
    mw
}
