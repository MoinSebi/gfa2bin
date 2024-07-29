use crate::core::core::MatrixWrapper;
use crate::core::helper::{merge_u32_to_u64, Feature};
use crate::window::window_main::getbv;
use bitvec::bitvec;
use bitvec::order::Lsb0;
use std::fmt::format;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};

use crate::block::block_main::block_wrapper;
use bitvec::prelude::BitVec;
use clap::ArgMatches;
use gfa_reader::{Gfa, Opt, Pansn, Segment};
use hashbrown::HashMap;
use log::info;
use rayon::prelude::*;

/// Subpath main function
///
/// Extract the subpath from a graph for each node
pub fn subpath_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // Read the arguments from the command line
    let graph_file = matches.value_of("gfa").unwrap();
    let output_prefix = matches.value_of("output").unwrap();
    let window: usize = matches.value_of("length").unwrap().parse().unwrap();
    let threads: usize = matches.value_of("threads").unwrap().parse().unwrap();
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
    subpath_wrapper(
        &wrapper,
        &graph,
        window,
        a,
        matches.is_present("blocks"),
        output_prefix,
        threads,
    )?;
    write_dummy_fam(&wrapper, &format!("{}.fam", output_prefix))?;
    Ok(())
}

/// Node to position index (hashmap)
///
/// For each path: Iterate over all nodes and store the index of each node
///
/// For each path:
///     node -> Vec<index>
///
/// Comment: Saved genome_id, haplotype_id, path_id, total_haplo, node2index
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
    blocks: bool,
    out_prefix: &str,
    threads: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    // Sample size
    let sample_size = graph2.genomes.len();
    let is_diploid = diploid_or_not(graph2);
    let segment_id = graph.segments.iter().map(|x| x.id).collect::<Vec<u32>>();
    segment_id
        .par_chunks(segment_id.len() / threads + 1)
        .enumerate()
        .for_each(|(i, chunks)| {
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

            for node_id in chunks {
                // Result vec
                let mut result_vec = get_traversals(*node_id, &node2index_hm, graph2, window);

                // Sort the traversals by genome_id, then haplotype_id
                result_vec.sort();

                // Checker
                if result_vec.is_empty() {
                    continue;
                }

                // !Thiis mmight be wring
                let vec_bitvec = traversal2bitvec(result_vec, sample_size, &is_diploid, &mut block);
                for x in 0..vec_bitvec.len() {
                    writeln!(
                        file_bim,
                        "{}\t{}\t{}",
                        "graph",
                        node_id.to_string() + &window.to_string() + &x.to_string(),
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
    concatenate_files_and_cleanup(
        &filenames1
            .iter()
            .map(|a1| format!("{}.bim", a1))
            .collect::<Vec<String>>(),
        format!("{}.bim", out_prefix),
        &Vec::new(),
    )?;
    concatenate_files_and_cleanup(
        &filenames1
            .iter()
            .map(|a1| format!("{}.bed", a1))
            .collect::<Vec<String>>(),
        format!("{}.bed", out_prefix),
        &vec![108, 27, 1],
    )?;
    if blocks {
        concatenate_files_and_cleanup(
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

pub fn make_filename(output_prefix: &str, num: usize) -> Vec<String> {
    (0..num)
        .into_iter()
        .map(|x| format!("{}_{}", output_prefix, x.to_string()))
        .collect()
}

pub fn get_traversals<'a>(
    node_id: u32,
    node2index_hm: &Vec<(usize, usize, usize, usize, HashMap<u32, Vec<usize>>)>,
    graph: &'a Pansn<u32, (), ()>,
    window: usize,
) -> Vec<(usize, usize, &'a [u32])> {
    let mut result_vec = Vec::new();
    for (genome_id, haplo_id, path_id, _total_haplo, node2index) in node2index_hm.iter() {
        // For each path, get the max number of nodes
        let max_index = graph.genomes[*genome_id].haplotypes[*haplo_id].paths[*path_id]
            .nodes
            .len();

        // If node is in path
        if node2index.contains_key(&node_id) {
            // Iterate over all occurrences
            for z in node2index.get(&node_id).unwrap() {
                // If window is within the path (does not exceed)
                // Push to vector with slice
                if window <= *z && *z + window < max_index + 1 {
                    result_vec.push((
                        *genome_id,
                        *haplo_id,
                        &graph.genomes[*genome_id].haplotypes[*haplo_id].paths[*path_id].nodes
                            [*z - window..*z + window],
                    ));
                }
            }
        }
    }
    result_vec
}

/// Convert collection of traversals to collection of samples with similar traversals
///
/// Genome_id, haplotype_id, path_id, &[u32]
///
/// Check if  &[u32] is the same
pub fn traversal2bitvec(
    traversals: Vec<(usize, usize, &[u32])>,
    number_of_samples: usize,
    ppl: &Vec<bool>,
    blocks: &mut Option<BufWriter<File>>,
) -> Vec<BitVec<u8>> {
    let mut sample_list: Vec<Vec<[usize; 2]>> = group_traversal(traversals);

    // If you want blocks written into extra
    if let Some(bufw) = blocks {
        writeln!(bufw, "{}\t{:?}", 0, sample_list).unwrap();
    }
    get_bitvector(&mut sample_list, number_of_samples, ppl)
}

/// Group traversals with similar traversals
///
/// One group contains all samples with the same traversal
/// Collected into one vector
pub fn group_traversal(traversals: Vec<(usize, usize, &[u32])>) -> Vec<Vec<[usize; 2]>> {
    let mut sample_list: Vec<Vec<[usize; 2]>> = Vec::new();
    let mut previous = traversals[0].2;
    sample_list.push(vec![[traversals[0].0, traversals[0].1]]);
    for traversal in traversals.iter().skip(1) {
        if traversal.2 != previous {
            previous = traversal.2;
            sample_list.push(vec![[traversal.0, traversal.1]]);
        } else {
            sample_list
                .last_mut()
                .unwrap()
                .push([traversal.0, traversal.1]);
        }
    }
    sample_list
}

/// Each group is one genotype.
///
/// Iterate over one genotype and setup a bitvector
///
pub fn get_bitvector(
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
use std::io::{self, BufReader, Read};
use std::path::Path;
use crate::core::bfile::write_dummy_fam;

pub fn concatenate_files_and_cleanup<P: AsRef<Path>>(
    input_files: &[P],
    output_file: P,
    buffer: &Vec<u8>,
) -> io::Result<()> {
    // Open the output file for writing
    let output = File::create(output_file)?;
    let mut writer = BufWriter::new(output);
    writer.write_all(buffer)?;
    // Process each input file
    for input_path in input_files {
        // Open the current input file for reading
        let input = File::open(input_path)?;
        let mut reader = BufReader::new(input);

        // Buffer to read data
        let mut buffer = [0; 8192]; // Adjust buffer size as needed

        // Read and write data in chunks
        loop {
            let bytes_read = reader.read(&mut buffer)?;
            if bytes_read == 0 {
                break; // End of file
            }
            writer.write_all(&buffer[..bytes_read])?;
        }

        // Delete the input file after processing
        std::fs::remove_file(input_path)?;
    }

    // Ensure all data is written to the output file
    writer.flush()?;

    Ok(())
}

pub fn diploid_or_not(gr: &Pansn<u32, (), ()>) -> Vec<bool> {
    let mut diploid_vec = Vec::new();
    for sample in gr.genomes.iter() {
        if sample.haplotypes.len() > 1 {
            diploid_vec.push(true);
        } else {
            diploid_vec.push(false);
        }
    }
    diploid_vec
}

