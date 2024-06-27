use crate::core::helper::{Feature};
use crate::r#mod::input_data::FileData;
use clap::ArgMatches;
use gfa_reader::{Gfa};



use crate::core::core::MatrixWrapper;







pub fn find_main(matches: &ArgMatches) {
    let graph_file = matches.value_of("gfa").unwrap();
    let feature_file = matches.value_of("features").unwrap();
    let _output = matches.value_of("output").unwrap();
    let _length = matches
        .value_of("length")
        .unwrap()
        .parse::<usize>()
        .unwrap();

    let data = FileData::from_file(feature_file);
    let feature = data.feature;

    if feature == Feature::MWindow {
        let mut mw = MatrixWrapper::new();
        mw.bfile_wrapper(feature_file);
        //find_matrix(data, mw, length);
    } else {
        let _graph: Gfa<u32, (), ()> = Gfa::parse_gfa_file(graph_file);
    }
}

// pub fn find_wrapper(
//     feature: Feature,
//     gfa: NCGfa<()>,
//     file_data: FileData,
//     length: usize,
//     output: &str,
// ) {
//     if feature == Feature::Node || feature == Feature::Edge || feature == Feature::DirNode {
//         find_easy(&gfa, file_data, feature, length);
//     } else if feature == Feature::Alignment {
//         find_alignment(&gfa, file_data, feature, length);
//     } else {
//         println!("not implemented")
//     }
// }

// pub fn find_easy(graph: &NCGfa<()>, input: FileData, feature: Feature, length: usize) {
//     let node_size = node_size(&graph);
//
//     //
//     let mut position_nodesize = Vec::new();
//     let mut vec_res_u64 = Vec::new();
//
//     for path in graph.paths.iter() {
//         let mut vec_u64 = Vec::new();
//         let mut index = Vec::new();
//         let mut pos = 0;
//         for i in 0..path.nodes.len() - 1 {
//             index.push([pos, node_size[path.nodes[i] as usize - 1]]);
//             pos += node_size[path.nodes[i] as usize - 1];
//             let v1 = path.nodes[i];
//             let v2 = path.dir[i];
//             let v3 = path.nodes[i + 1];
//             let v4 = path.dir[i + 1];
//
//             if feature == Feature::Node {
//                 vec_u64.push(v1 as u64);
//             } else if feature == Feature::DirNode {
//                 vec_u64.push(v1 as u64 * 2 + v2 as u64);
//             } else if feature == Feature::Edge {
//                 let u1 = v1 * 2 + v2 as u32;
//                 let u2 = v3 * 2 + v4 as u32;
//                 vec_u64.push(merge_u32_to_u64(u1, u2));
//             }
//         }
//         vec_res_u64.push(vec_u64);
//         position_nodesize.push(index)
//     }
//
//     let file = File::create("output").unwrap();
//     let mut writer = std::io::BufWriter::new(file);
//     let data_hs = input.data.iter().collect::<HashSet<&u64>>();
//     for (i, x) in vec_res_u64.iter().enumerate() {
//         for (i2, y) in x.iter().enumerate() {
//             if data_hs.contains(y) {
//                 writeln!(
//                     writer,
//                     "{}\t{}\t{}\tID:{};NS:{};NB:{}",
//                     graph.paths[i].name,
//                     max(0, position_nodesize[i][i2][0] as i128 - length),
//                     position_nodesize[i][i2][0] as i128
//                         + position_nodesize[i][i2][1] as i128
//                         + 1
//                         + length,
//                     to_string1(*y, &feature),
//                     position_nodesize[i][i2][1],
//                     position_nodesize[i][i2][0],
//                 )
//                 .expect("Error writing to file")
//             }
//         }
//     }
// // }
//
// pub fn find_alignment(graph: &NCGfa<()>, input: FileData, feature: Feature, length: usize) {
//     let mut ff = FileData::new();
//     ff.data = input
//         .data
//         .iter()
//         .map(|x| split_u64_to_u32s(*x).0 as u64)
//         .collect();
//     ff.feature = Feature::Node;
//     find_easy(graph, ff, feature, length);
// }
//
// pub fn find_matrix(input: FileData, mw: MatrixWrapper, window: usize) {
//     let num_path = mw.matrix_bit[0].len() / 2;
//
//     let mut j = 0;
//     for x in window..mw.matrix_bit.len() - window {
//         if split_u64_to_u32s(input.data[j]).0 as u64 == mw.geno_names[x] {
//             j += 1;
//             let ff = split_u64_to_u32s(input.data[j]).1;
//             let mut bv2 = Vec::new();
//             for y in 0..num_path {
//                 let mut bv: Vec<[bool; 2]> = Vec::new();
//
//                 for z in x - window..x + window {
//                     bv.push([mw.matrix_bit[z][y * 2], mw.matrix_bit[z][y * 2 + 1]]);
//                 }
//                 bv2.push((y, bv));
//             }
//             // sort by the bitvec
//             bv2.sort_by(|a, b| a.1.cmp(&b.1));
//             let f = &get_index(&bv2)[ff as usize];
//         }
//     }
// }
//
// pub fn write_matrix(
//     mw: MatrixWrapper,
//     output: &str,
//     win: usize,
//     start: usize,
//     aa: Vec<(usize, Vec<[bool; 2]>)>,
// ) {
//     let file = File::create(output).unwrap();
//     let mut writer = std::io::BufWriter::new(file);
//
//     write!(writer, "{}", mw.geno_names[start - win, start + win]).expect("Error writing to file");
//     for x in aa.iter() {
//         write!(writer, "{}", mw.sample_names[x.0]).expect("Error writing to file");
//         write!(
//             writer,
//             "{}",
//             x.1.iter().map(|x| matric_writer(x)).collect::<String>()
//         )
//         .expect("Error writing to file");
//     }
// }
//
// pub fn matric_writer(a: &[bool; 2]) -> String {
//     return format!("{}{}", a[0] as u8, a[1] as u8);
// }

// pub fn find_subpath(
//     graph2: &Pansn<NCPath>,
//     nodes: Vec<u32>,
//     window: usize,
//     vv: Vec<(usize, usize, usize, HashMap<u32, Vec<usize>>)>,
// ) {
//     let sample_size = graph2.genomes.len();
//     for x1 in nodes.iter() {
//         let mut vecc = Vec::new();
//         for (genome_id, haplotype_id, path_id, node2index) in vv.iter() {
//             let ll = graph2.genomes[*genome_id].haplotypes[*haplotype_id].paths[*path_id]
//                 .nodes
//                 .len();
//             if node2index.contains_key(&x1) {
//                 for z in node2index.get(&x1).unwrap() {
//                     if window <= *z && *z + window < ll + 1 {
//                         vecc.push((
//                             graph2.genomes[*genome_id].haplotypes[*haplotype_id].paths[*path_id]
//                                 .nodes[*z - window..*z + window],
//                             *genome_id,
//                             *haplotype_id,
//                             *path_id,
//                             *z,
//                             *window,
//                         ));
//                     }
//                 }
//             }
//         }
//         vecc.sort();
//         if vecc.is_empty() {
//             continue;
//         }
//     }
// }
//
// pub fn write_subpath(
//     aa: &Vec<(Vec<&[u32]>, usize, usize, usize, usize, usize)>,
//     ll: usize,
//     output: &str,
//     pp: &Pansn<NCPath>,
// ) {
//     let mut ab = aa.first().unwrap();
//     let file = File::create(output).unwrap();
//     let mut writer = std::io::BufWriter::new(file);
//     let mut c = 0;
//
//     for y in aa.iter() {
//         if y.0 != ab.0 {
//             c += 1;
//             ab = y;
//         }
//
//         if c == ll {
//             write!(
//                 writer,
//                 "{}",
//                 pp.genomes[ab.1].haplotypes[ab.2].paths[ab.3].name
//             )
//             .expect("Error writing to file");
//             write!(writer, "{}", 1).expect("Error writing to file");
//         }
//     }
// }
