use crate::core::core::MatrixWrapper;
use crate::core::helper::{merge_u32_to_u64, Feature};
use crate::subpath::subpath_main::{diploid_or_not, make_filename, traversal2bitvec};
use clap::ArgMatches;
use gfa_reader::{Gfa, Pansn};
use rayon::prelude::*;

use log::info;
use std::collections::HashSet;
use std::fmt::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::{fmt, io};
use crate::core::bfile::write_dummy_fam;

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
    let threads: usize = matches.value_of("threads").unwrap().parse().unwrap();

    info!("Block subcommand");
    info!("Graph file: {}", graph_file);
    info!("Separator: {}", sep);
    info!("Window size: {}", window);
    info!("Distance: {}", cutoff_distance);
    info!("Step size: {}", step_size);
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
    let aa = wrapper_blocks(&wrapper, b, a, cutoff_distance, true, output_prefix, threads)?;
    write_dummy_fam(&wrapper, &format!("{}.fam", output_prefix))?;
    Ok(())
}

pub fn block_wrapper(
    graph: &Gfa<u32, (), ()>,
    step: usize,
    window_size: usize,
    sequence: Option<&str>,
    sequence_window: Option<&str>,
) -> Result<Vec<[u32; 2]>, Box<dyn std::error::Error>> {
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
    for x in (0..(graph.segments.len() - wsize)).step_by(step) {
        gg.push([graph.segments[x].id, graph.segments[x + wsize].id]);
    }
    gg
}

pub fn block_seq(graph: &Gfa<u32, (), ()>, size: usize, step_size: usize) -> Vec<[u32; 2]> {
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
    blocks: bool,
    out_prefix: &str,
    threads: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let sample_size = graph2.genomes.len();
    let is_diploid = diploid_or_not(graph2);

    block
        .par_chunks(block.len() / threads + 1)
        .enumerate()
        .for_each(|(i, chunk)| {
            let mut file_bed = BufWriter::new(
                File::create(format!("{}_{}.bed", out_prefix, i))
                    .expect("Not able to create bed file"),
            );
            let mut file_bim = BufWriter::new(
                File::create(format!("{}_{}.bim", out_prefix, i))
                    .expect("Not able to create bim file"),
            );
            let mut block = None;
            if blocks {
                block = Some(BufWriter::new(
                    File::create(format!("{}_{}.block", out_prefix, i))
                        .expect("Not able to create block file"),
                ));
            }

            for block_border in chunk.iter() {
                println!("Block: {:?}", block_border);
                let block_hashset = (block_border[0]..block_border[1]).collect::<HashSet<u32>>();
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
                                    &x.nodes[block_array[0]..block_array[1]],
                                ));
                            }
                        }
                    }
                }
                all_blocks.sort();
                if all_blocks.is_empty() {
                    continue;
                }
                let vec_bitvec = traversal2bitvec(all_blocks, sample_size, &is_diploid, &mut block);
                for x in 0..vec_bitvec.len() {
                    writeln!(
                        file_bim,
                        "{}\t{}\t{}",
                        "graph",
                        block_border[0].to_string() + "_" + &block_border[1].to_string() + &x.to_string(),
                        x
                    )
                    .unwrap();
                    let buff = vec_bitvec[x].as_raw_slice();
                    file_bed.write_all(&buff).expect("Not able to write ")
                }
            }
        });
    info!("Concatenating files");
    let filenames1 = make_filename(out_prefix, threads);
    crate::subpath::subpath_main::concatenate_files_and_cleanup(
        &filenames1
            .iter()
            .map(|a1| format!("{}.bim", a1))
            .collect::<Vec<String>>(),
        format!("{}.bim", out_prefix),
        &Vec::new(),
    )?;
    crate::subpath::subpath_main::concatenate_files_and_cleanup(
        &filenames1
            .iter()
            .map(|a1| format!("{}.bed", a1))
            .collect::<Vec<String>>(),
        format!("{}.bed", out_prefix),
        &vec![108, 27, 1],
    )?;
    if blocks {
        crate::subpath::subpath_main::concatenate_files_and_cleanup(
            &filenames1
                .iter()
                .map(|a1| format!("{}.block", a1))
                .collect::<Vec<String>>(),
            format!("{}.block", out_prefix),
            &Vec::new(),
        )?;
    }
    info!("Done");
    Ok(())
}
