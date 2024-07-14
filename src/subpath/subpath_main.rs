use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use crate::core::core::MatrixWrapper;
use crate::core::helper::{merge_u32_to_u64, Feature};
use crate::window::window_main::getbv;
use bitvec::bitvec;
use bitvec::order::Lsb0;

use bitvec::prelude::BitVec;
use clap::ArgMatches;
use gfa_reader::{Gfa, Opt, Pansn, Segment};
use hashbrown::HashMap;
use log::info;

/// Subpath main function
///
/// Extract the subpath from a graph for each node
pub fn subpath_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // Read the arguments from the command line
    let graph_file = matches.value_of("gfa").unwrap();
    let output_prefix = matches.value_of("output").unwrap();
    let window: usize = matches.value_of("length").unwrap().parse().unwrap();
    let split = matches
        .value_of("split")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();

    let mut block = None;
    if let Some(blocks_path) = matches.value_of("blocks") {
        let file = File::create(blocks_path).unwrap();
        block = Some(BufWriter::new(file));
    }

    /// Check the arguments
    info!("Subpath subcommand");
    info!("Graph file: {}", graph_file);
    info!("Output prefix: {}", output_prefix);
    info!("Window length: {}", window);

    info!("Reading graph file");
    let mut graph: Gfa<u32, (), ()> = Gfa::parse_gfa_file(graph_file);
    graph.walk_to_path("#");
    let wrapper: Pansn<u32, (), ()> = Pansn::from_graph(&graph.paths, "#");

    info!("Indexing graph");
    let a = gfa_index(&wrapper);

    info!("Extracting subpath");
    let mw = subpath_wrapper(&wrapper, &graph, window, a, &mut block);

    info!("Number of nodes: {}", graph.segments.len());
    info!("Number of subpath: {}", mw.matrix_bit.len());
    info!("Number of subpath: {}", mw.geno_names.len());

    info!("Writing output");

    Ok(())
}

/// #Index the graph
///
/// For each path:
///     node -> Vec<index>
pub fn gfa_index(
    graph: &Pansn<u32, (), ()>,
) -> Vec<(usize, usize, usize, usize, HashMap<u32, Vec<usize>>)> {
    // Index (genome_id, haplotype_id, path_id, {node -> Vec<index>})
    let mut index = Vec::new();
    for (genome_id, path) in graph.genomes.iter().enumerate() {
        for (haplo_id, x) in path.haplotypes.iter().enumerate() {
            for (path_id, p) in x.paths.iter().enumerate() {
                let mut node2index: HashMap<u32, Vec<usize>> = HashMap::new();

                for (ind, node) in p.nodes.iter().enumerate() {
                    node2index
                        .entry(node.clone())
                        .or_insert_with(Vec::new)
                        .push(ind);
                }
                index.push((
                    genome_id,
                    haplo_id,
                    path_id,
                    path.haplotypes.len(),
                    node2index,
                ));
            }
        }
    }
    index
}

/// Extract all subpath for each node
///
pub fn subpath_wrapper(
    graph2: &Pansn<u32, (), ()>,
    graph: &Gfa<u32, (), ()>,
    window: usize,
    node2index_hm: Vec<(usize, usize, usize, usize, HashMap<u32, Vec<usize>>)>,
    blocks: &mut Option<BufWriter<File>>,
) -> MatrixWrapper {
    // Initialize the matrix wrapper
    let mut mw = MatrixWrapper::new();

    // Sample size
    let sample_size = graph2.genomes.len();
    let is_diploid = diploid_or_not(graph2);
    // Iterate over each node
    let a = graph.segments.iter().map(|x| x.id).collect::<Vec<u32>>();

    let chunks = a.chunks(a.len() / 8);

    for node_id in graph.segments.iter() {
        // Result vec
        let mut result_vec = get_t(&node_id, &node2index_hm, graph2, window);


        // Sort by first then second thing
        result_vec.sort();
        if result_vec.is_empty() {
            continue;
        }

        // !Thiis mmight be wring
        let o = traversal2samples(result_vec, sample_size, &is_diploid, blocks);

        mw.matrix_bit.extend(o);
    }
    println!("{:?}", mw.matrix_bit.len());
    mw.feature = Feature::PWindow;
    mw.window_size = window;
    mw.shape = (mw.matrix_bit.len(), mw.matrix_bit[0].len());
    mw.sample_names = graph2.genomes.iter().map(|x| x.name.clone()).collect();

    mw
}

pub fn get_t<'a>(node_id: &Segment<u32, ()>, node2index_hm: &Vec<(usize, usize, usize, usize, HashMap<u32, Vec<usize>>)>,   graph2: &'a  Pansn<u32, (), ()>,
             window: usize) -> Vec<(usize, usize, &'a [u32])> {
    let mut result_vec = Vec::new();
    for (genome_id, haplotype_id, path_id, total_haplo, node2index) in node2index_hm.iter() {
        // Max index
        let max_index = graph2.genomes[*genome_id].haplotypes[*haplotype_id].paths[*path_id]
            .nodes
            .len();

        //
        if node2index.contains_key(&node_id.id) {
            // Iterate of there are more than one occurence of the node
            for z in node2index.get(&node_id.id).unwrap() {
                if window <= *z && *z + window < max_index + 1 {
                    result_vec.push((
                        *genome_id,
                        *haplotype_id,
                        &graph2.genomes[*genome_id].haplotypes[*haplotype_id].paths[*path_id]
                            .nodes[*z - window..*z + window],
                    ));
                }
            }
        }
    }
    result_vec
}

pub fn diploid_or_not(gr: &Pansn<u32, (), ()>) -> Vec<bool> {
    let mut ii = Vec::new();
    for x in gr.genomes.iter() {
        if x.haplotypes.len() > 1 {
            ii.push(true);
        } else {
            ii.push(false);
        }
    }
    ii
}

/// Convert collection of traversals to collection of samples with similar traversals
///
/// Genome_id, haplotype_id, path_id, &[u32]
///
/// Check if  &[u32] is the same
pub fn traversal2samples(
    ii: Vec<(usize, usize, &[u32])>,
    number_of_samples: usize,
    ppl: &Vec<bool>,
    blocks: &mut Option<BufWriter<File>>
) -> Vec<BitVec<u8>> {
    let mut pp: Vec<Vec<[usize; 2]>> = Vec::new();
    let mut last = ii[0].2;
    pp.push(vec![[ii[0].0, ii[0].1]]);
    for x in ii.iter().skip(1) {
        if x.2 != last {
            last = x.2;
            pp.push(vec![[x.0, x.1]]);
        } else {
            pp.last_mut().unwrap().push([x.0, x.1]);
        }
    }
    if let Some(b) = blocks {
        writeln!(b, "{}\t{:?}", 0, pp).unwrap();
    }
    getbv2(&mut pp, number_of_samples, ppl)
}

/// Create a bitvector
pub fn getbv2(
    present_sample_collection: &mut Vec<Vec<[usize; 2]>>,
    len: usize,
    is_diploid: &Vec<bool>,
) -> Vec<BitVec<u8>> {
    let mut bitvec_collection = Vec::new();
    for samples in present_sample_collection.iter_mut() {
        let mut bitvec_tmp: BitVec<u8, Lsb0> = BitVec::<u8, Lsb0>::repeat(false, len * 2 + 1);
        let mut index = 0;
        while index < samples.len() - 1 {
            if samples[index][0] == samples[index + 1][0] {
                bitvec_tmp.set(samples[index][0] * 2, true);
                bitvec_tmp.set(samples[index][0] * 2 + 1, true);
                index += 2;
            } else {
               if !is_diploid[samples[index][0]] {
                    bitvec_tmp.set(samples[index][0] * 2, true);
                    bitvec_tmp.set(samples[index][0] * 2 + 1, true);
                    index += 1;
                } else {
                    bitvec_tmp.set(samples[index][0] * 2 + 1, true);
                    index += 1;
                }
            }

        }
        if index == samples.len() - 1 {
            bitvec_tmp.set(samples[index][0] * 2 + 1, true);
        }
        bitvec_collection.push(bitvec_tmp);

    }
    bitvec_collection
}
