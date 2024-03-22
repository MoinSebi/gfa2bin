use clap::ArgMatches;
use gfa_reader::{NCGfa, NCPath, Pansn};
use hashbrown::HashMap;
use bitvec::order::Lsb0;
use bitvec::prelude::BitVec;
use log::info;
use crate::core::core::MatrixWrapper;
use crate::core::helper::{Feature, merge_u32_to_u64};
use crate::window::window_main::getbv;


/// Subpath main function
///
/// Extract the subpath from a graph for each node
pub fn subpath_main(matches: &ArgMatches) {
    let graph_file = matches.value_of("gfa").unwrap();
    let output_prefix = matches.value_of("output").unwrap();
    let window: usize = matches.value_of("length").unwrap().parse().unwrap();
    let split = matches
        .value_of("split")
        .unwrap_or("1")
        .parse::<usize>()
        .unwrap();

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
    let mw = wrapper_matrix(&wrapper, &graph, window, a);

    info!("Number of nodes: {}", graph.nodes.len());
    info!("Number of subpath: {}", mw.matrix_bit.len());
    let feature = Feature::PWindow;

    info!("Writing output");
    let chunk_size = (mw.matrix_bit.len() / split) + 1;
    let chunks = mw.matrix_bit.chunks(chunk_size);

    let len = chunks.len();
    for (index, _y) in chunks.enumerate() {
        mw.write_fam(index, output_prefix, feature, len);
        mw.write_bed(index, output_prefix, feature, len);
        mw.write_bim(index, output_prefix, &feature, len);
    }
}

/// Index the graph
///
/// For each path:
///     node -> index
pub fn gfa_index(graph: &Pansn<NCPath>) -> Vec<(usize, usize, usize, HashMap<u32, Vec<usize>>)> {
    let mut vv = Vec::new();
    for (genome_id, path) in graph.genomes.iter().enumerate() {
        for (haplotype_id, x) in path.haplotypes.iter().enumerate() {
            for (path_id, p) in x.paths.iter().enumerate() {
                let mut node2index: HashMap<u32, Vec<usize>> = HashMap::new();

                for (i, y) in p.nodes.iter().enumerate() {
                    if node2index.contains_key(y) {
                        node2index.get_mut(y).unwrap().push(i);
                    } else {
                        node2index.insert(*y, vec![i]);
                    }
                }
                vv.push((genome_id, haplotype_id, path_id, node2index));
            }
        }
    }

    return vv;
}

/// Extract all subpath for each node
pub fn wrapper_matrix(
    graph2: &Pansn<NCPath>,
    graph: &NCGfa<()>,
    window: usize,
    vv: Vec<(usize, usize, usize, HashMap<u32, Vec<usize>>)>,
) -> MatrixWrapper{
    let mut mw = MatrixWrapper::new();
    let difference = graph2.genomes.len();
    for x1 in graph.nodes.iter() {
        let mut vecc = Vec::new();
        for (genome_id, haplotype_id, path_id, node2index) in vv.iter() {
            let ll = graph2.genomes[*genome_id].haplotypes[*haplotype_id].paths[*path_id].nodes.len();
            if node2index.contains_key(&x1.id) {
                for z in node2index.get(&x1.id).unwrap() {
                    if window <= *z &&  *z + window < ll + 1{
                        vecc.push((
                            genome_id,
                            haplotype_id,
                            path_id,
                            &graph2.genomes[*genome_id].haplotypes[*haplotype_id].paths[*path_id].nodes
                                [*z - window..
                                *z + window],
                        ));
                    }
                }
            }
        }
        vecc.sort();
        if vecc.is_empty() {
            continue;
        }
        mw.geno_names.extend((0..vecc.len()).map(|_| merge_u32_to_u64(x1.id, vecc.len() as u32)));

        mw.matrix_bit.extend(function1(vecc, difference));
    }
    mw.feature = Feature::PWindow;
    mw.shape = (mw.matrix_bit.len(), mw.matrix_bit[0].len());
    mw.fam_entries = graph2.genomes.iter().map(|x| x.name.clone()).collect();


    return mw
}

pub fn function1(ii: Vec<(&usize, &usize, &usize, &[u32])>, difference: usize) -> Vec<BitVec<u8>>{
    let mut pp = Vec::new();
    let mut last = ii[0].3;
    pp.push(vec![*ii[0].0]);
    for x in ii.iter().skip(1) {
        if x.3 != last {
            last = x.3;
            pp.push(vec![*x.0]);
        } else {
            pp.last_mut().unwrap().push(*x.0);
        }
    }
    getbv(&pp, difference)
}

