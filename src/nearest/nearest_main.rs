use clap::ArgMatches;
use gfa_reader::Gfa;
use hashbrown::{HashMap, HashSet};
use log::info;
use std::fs::File;

use std::i64;
use std::io::{self, BufWriter};
use std::io::{BufRead, BufReader, Write};
use std::str::FromStr;

/// Nearest node main function
///
/// Find the closest reference node for each node in the graph
pub fn nearest_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    info!("Running 'gfa2bin nearest'");

    // Graph file
    let graph_file = matches.value_of("gfa").unwrap();

    // Output file
    let output_file = matches.value_of("output").unwrap();

    let threads = matches.value_of("threads").unwrap();

    info!("GFA file: {}", graph_file);
    info!(
        "Node file: {}",
        if matches.is_present("nodes") {
            matches.value_of("nodes").unwrap()
        } else {
            "All nodes"
        }
    );
    info!(
        "Reference file: {}",
        if matches.is_present("references") {
            matches.value_of("references").unwrap()
        } else {
            "None"
        }
    );
    info!(
        "Prefix: {}",
        if matches.is_present("prefix") {
            matches.value_of("prefix").unwrap()
        } else {
            "None"
        }
    );
    info!("Threads: {}", threads);
    info!("Output file: {}\n", output_file);

    info!("Read GFA file");
    let mut graph = Gfa::parse_gfa_file_multi(graph_file, threads.parse().unwrap());

    info!("Convert walks to path with '#' separator");
    graph.walk_to_path("#");

    // Requested nodes
    let requested_nodes: Vec<u32>;
    if matches.is_present("nodes") {
        info!(
            "Reading requested nodes from {}",
            matches.value_of("nodes").unwrap()
        );
        requested_nodes = read_input(matches.value_of("nodes").unwrap()).unwrap();
    } else {
        info!("No nodes provided, all nodes will be considered");
        requested_nodes = graph.segments.iter().map(|x| x.id).collect();
    }

    // Which path are "reference" paths
    let mut ref_list = Vec::new();
    if matches.is_present("references") {
        ref_list = read_input(matches.value_of("references").unwrap())?;
    } else if matches.is_present("prefix") {
        ref_list = by_prefix(matches.value_of("prefix").unwrap(), &graph)?;
    } else {
        panic!("You need to provide either a reference list or a prefix")
    }

    info!("Finding closest reference node for each node");
    let closest_node_vec = read_nodes(
        &graph,
        &ref_list,
        &requested_nodes.iter().cloned().collect::<HashSet<u32>>(),
    );

    info!("Writing output to {}", output_file);
    write_file(closest_node_vec, output_file, &graph, &ref_list).unwrap();

    Ok(())
}

/// # Extract reference paths by prefix
///
/// - .startswith(prefix) on all path
pub fn by_prefix(
    prefix: &str,
    graph: &Gfa<u32, (), ()>,
) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let reference_paths = graph
        .paths
        .iter()
        .map(|x| x.name.clone())
        .filter(|x| x.starts_with(prefix))
        .collect::<Vec<_>>();
    if reference_paths.is_empty() {
        Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "No reference paths found",
        )))
    } else {
        info!("{:?}", reference_paths);
        Ok(reference_paths)
    }
}

/// Read input by line
///
/// Read input by line and return a vector of the type T
pub fn read_input<T: FromStr>(input: &str) -> Result<Vec<T>, io::Error> {
    let file = File::open(input).unwrap_or_else(|_| panic!("Can not open file: {}", input));
    let reader = BufReader::new(file);

    let mut lines_vec: Vec<T> = Vec::new();

    for line in reader.lines() {
        let line = line?; // This unwraps Result<String, io::Error> or returns an error if it occurs
        match line.parse::<T>() {
            Ok(value) => lines_vec.push(value),
            Err(_) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Failed to parse node (line)",
                ))
            }
        }
    }
    Ok(lines_vec)
}

/// Get a (node, position) vector for a specific path
pub fn pos(graph: &Gfa<u32, (), ()>, name: String) -> Vec<(u32, usize)> {
    for path in graph.paths.iter() {
        if path.name == name {
            let mut bb = Vec::new();
            let mut pos1: usize = 0;

            for node in path.nodes.iter() {
                bb.push((*node, pos1));
                pos1 += graph.get_sequence_by_id(node).len();
            }
            bb.sort_by(|a, b| a.0.cmp(&b.0));
            return bb;
        }
    }
    Vec::new()
}

/// # Collect reference nodes
///
/// - All nodes in reference paths
/// - Collection is a hashset for better testing
pub fn get_ref_nodes(graph: &Gfa<u32, (), ()>, names: &Vec<String>) -> HashSet<u32> {
    let mut nodes_hashset = HashSet::new();
    for path in graph.paths.iter() {
        if names.contains(&path.name) {
            nodes_hashset.extend(&path.nodes)
        }
    }
    nodes_hashset
}

/// # Initialize a hashmap of node -> [ref_node, distance]
///
/// This can only be run for one reference at once
pub fn init_hm(graph: &Gfa<u32, (), ()>, checked_nodes: &HashSet<u32>) -> HashMap<u32, [i64; 2]> {
    let mut result_hm = HashMap::with_capacity(graph.segments.len());
    if checked_nodes.is_empty() {
        for x in graph.segments.iter() {
            result_hm.insert(x.id, [x.id as i64, i64::MAX]);
        }
    } else {
        for x in checked_nodes.iter() {
            result_hm.insert(*x, [*x as i64, i64::MAX]);
        }
    }
    result_hm
}

/// For each node in the graph, find the closest reference node and the distance
///
/// Return a vector [(node, ref_node, distance)]
pub fn read_nodes(
    graph: &Gfa<u32, (), ()>,
    names: &Vec<String>,
    checked_nodes: &HashSet<u32>,
) -> Vec<(u32, i64, i64)> {
    let mut distance = 0;
    let mut reference_node = 0;
    let nodes_hs = get_ref_nodes(graph, names);
    let mut result_hm = init_hm(graph, checked_nodes);
    for path in graph.paths.iter() {
        if !names.contains(&path.name) {
            for node in path.nodes.iter() {
                if nodes_hs.contains(node) {
                    distance = 0;
                    reference_node = *node;
                    *result_hm.get_mut(node).unwrap() = [*node as i64, -1]
                } else {
                    if checked_nodes.contains(node) && result_hm[node][1] > distance {
                        *result_hm.get_mut(node).unwrap() = [reference_node as i64, distance]
                    }
                    distance += graph.get_sequence_by_id(node).len() as i64;
                }
            }
            for node in path.nodes.iter().rev() {
                if nodes_hs.contains(node) {
                    distance = 0;
                    reference_node = *node;
                    *result_hm.get_mut(node).unwrap() = [*node as i64, -1]
                } else {
                    if checked_nodes.contains(node) && result_hm[node][1] > distance {
                        *result_hm.get_mut(node).unwrap() = [reference_node as i64, distance]
                    }
                    distance += graph.get_sequence_by_id(node).len() as i64;
                }
            }
        }
    }

    // for reference nodes
    for path in graph.paths.iter() {
        if names.contains(&path.name) {
            for node in path.nodes.iter() {
                *result_hm.get_mut(node).unwrap() = [*node as i64, -1]
            }
        }
    }

    // Write result
    let mut result_vec = result_hm
        .iter()
        .map(|(k, v)| (*k, v[0], v[1]))
        .collect::<Vec<_>>();
    result_vec.sort_by(|a, b| a.1.cmp(&b.1));
    result_vec
}

/// # Write 'nearest' output to file
///
/// Double index approach for fast processing
pub fn write_file(
    result: Vec<(u32, i64, i64)>,
    output: &str,
    graph: &Gfa<u32, (), ()>,
    names: &Vec<String>,
) -> Result<(), std::io::Error> {
    // Create file
    let file_out = File::create(output).expect("Unable to create file");
    let mut output_reader = BufWriter::new(file_out);
    writeln!(output_reader, "node\tref_node\tdistance\tposition\tpath")?;

    for reference_name in names.iter() {
        // Have multiple pos for a single node
        let pos1 = pos(graph, reference_name.to_string());
        let mut i = 0;
        let mut j = 0;

        while i < result.len() && j < pos1.len() {
            if result[i].1 == pos1[j].0 as i64 {
                let value = result[i].1;

                // Find the range of indices for the matching value in vec1
                let start_i = i;
                while i < result.len() && result[i].1 == value {
                    i += 1;
                }

                // Find the range of indices for the matching value in vec2
                let start_j = j;
                while j < pos1.len() && pos1[j].0 as i64 == value {
                    j += 1;
                }

                // Directly generate combinations of indices for the matching value
                for index1 in start_i..i {
                    for index2 in start_j..j {
                        writeln!(
                            output_reader,
                            "{}\t{}\t{}\t{}\t{}",
                            result[index1].0,
                            result[index1].1,
                            result[index1].2,
                            pos1[index2].1,
                            reference_name
                        )?;
                    }
                }
            } else if result[i].1 < pos1[j].0 as i64 {
                i += 1;
            } else {
                j += 1;
            }
        }
    }
    Ok(())
}
