use crate::core::core::MatrixWrapper;
use crate::core::helper::{merge_u32_to_u64, Feature};
use crate::window::window_main::getbv;
use bitvec::order::Lsb0;
use bitvec::prelude::BitVec;
use clap::ArgMatches;
use gfa_reader::{NCGfa, NCPath, Pansn};
use hashbrown::HashMap;
use log::info;

/// Subpath main function
///
/// Extract the subpath from a graph for each node
pub fn subpath_main(matches: &ArgMatches) {
    // Read the arguments from the command line
    let graph_file = matches.value_of("gfa").unwrap();
    let output_prefix = matches.value_of("output").unwrap();
    let window: usize = matches.value_of("length").unwrap().parse().unwrap();
    let split = matches
        .value_of("split")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();

    /// Check the arguments
    info!("Subpath subcommand");
    info!("Graph file: {}", graph_file);
    info!("Output prefix: {}", output_prefix);
    info!("Window length: {}", window);

    info!("Reading graph file");
    let mut graph: NCGfa<()> = NCGfa::new();
    graph.parse_gfa_file(graph_file, false);
    let wrapper: Pansn<NCPath> = Pansn::from_graph(&graph.paths, "#");

    info!("Indexing graph");
    let a = gfa_index(&wrapper);

    info!("Extracting subpath");
    let mw = subpath_wrapper(&wrapper, &graph, window, a);

    info!("Number of nodes: {}", graph.nodes.len());
    info!("Number of subpath: {}", mw.matrix_bit.len());

    info!("Writing output");
    mw.write_chunks(split, output_prefix, Feature::PWindow);
}

/// Index the graph
///
/// For each path:
///     node -> Vec<index>
pub fn gfa_index(graph: &Pansn<NCPath>) -> Vec<(usize, usize, usize, HashMap<u32, Vec<usize>>)> {
    let mut index = Vec::new();
    for (genome_id, path) in graph.genomes.iter().enumerate() {
        for (haplo_id, x) in path.haplotypes.iter().enumerate() {
            for (path_id, p) in x.paths.iter().enumerate() {
                let mut node2index: HashMap<u32, Vec<usize>> = HashMap::new();

                for (i, y) in p.nodes.iter().enumerate() {
                    if node2index.contains_key(y) {
                        node2index.get_mut(y).unwrap().push(i);
                    } else {
                        node2index.insert(*y, vec![i]);
                    }
                }
                index.push((genome_id, haplo_id, path_id, node2index));
            }
        }
    }

    return index;
}

/// Extract all subpath for each node
pub fn subpath_wrapper(
    graph2: &Pansn<NCPath>,
    graph: &NCGfa<()>,
    window: usize,
    vv: Vec<(usize, usize, usize, HashMap<u32, Vec<usize>>)>,
) -> MatrixWrapper {
    let mut mw = MatrixWrapper::new();
    let sample_size = graph2.genomes.len();
    for x1 in graph.nodes.iter() {
        let mut vecc = Vec::new();
        for (genome_id, haplotype_id, path_id, node2index) in vv.iter() {
            let ll = graph2.genomes[*genome_id].haplotypes[*haplotype_id].paths[*path_id]
                .nodes
                .len();
            if node2index.contains_key(&x1.id) {
                for z in node2index.get(&x1.id).unwrap() {
                    if window <= *z && *z + window < ll + 1 {
                        vecc.push((
                            *genome_id,
                            *haplotype_id,
                            *path_id,
                            &graph2.genomes[*genome_id].haplotypes[*haplotype_id].paths[*path_id]
                                .nodes[*z - window..*z + window],
                        ));
                    }
                }
            }
        }
        vecc.sort();
        if vecc.is_empty() {
            continue;
        }
        mw.geno_names
            .extend((0..vecc.len()).map(|a| merge_u32_to_u64(x1.id, a as u32)));

        mw.matrix_bit.extend(function1(vecc, sample_size));
    }
    mw.feature = Feature::PWindow;
    mw.window_size = window;
    mw.shape = (mw.matrix_bit.len(), mw.matrix_bit[0].len());
    mw.sample_names = graph2.genomes.iter().map(|x| x.name.clone()).collect();

    return mw;
}

pub fn function1(
    ii: Vec<(usize, usize, usize, &[u32])>,
    number_of_samples: usize,
) -> Vec<BitVec<u8>> {
    let mut pp = Vec::new();
    let mut last = ii[0].3;
    pp.push(vec![ii[0].0]);
    for x in ii.iter().skip(1) {
        if x.3 != last {
            last = x.3;
            pp.push(vec![x.0]);
        } else {
            pp.last_mut().unwrap().push(x.0);
        }
    }
    getbv(&pp, number_of_samples)
}
