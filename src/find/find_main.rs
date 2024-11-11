use crate::core::helper::merge_u32_to_u64;
use clap::ArgMatches;
use gfa_reader::Gfa;
use std::cmp::PartialEq;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};

/// Main function for find subcommand
pub fn find_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    let graph_file = matches.value_of("gfa").unwrap();
    let feature_file = matches.value_of("features").unwrap();
    let output = matches.value_of("output").unwrap();
    let _length = matches
        .value_of("length")
        .unwrap()
        .parse::<usize>()
        .unwrap();

    let a = determine_type(feature_file)?;
    find_easy(
        &Gfa::parse_gfa_file(graph_file),
        &a,
        read_file_lines(feature_file, &a)?,
        output,
    )?;
    Ok(())
}

#[derive(PartialEq, Eq, Hash)]
enum InputType {
    Segment,
    DirSegment,
    Link,
    Subgraph,
    Block,
}

pub fn to_string(dig: u64, input_type: InputType) -> String {
    if input_type == InputType::Segment {
        dig.to_string()
    } else if input_type == InputType::DirSegment {
        let dig1 = dig / 2;
        let dig2 = dig % 2;
        return format!("{}{}", dig1, if dig2 == 0 { "+" } else { "-" });
    } else if input_type == InputType::Link {
        let dig1 = dig / 2;
        let dig2 = dig % 2;
        let dig3 = dig1 / 2;
        let dig4 = dig1 % 2;
        return format!(
            "{}{}{}{}",
            dig3,
            if dig4 == 0 { "+" } else { "-" },
            dig2 / 2,
            if dig2 % 2 == 0 { "+" } else { "-" }
        );
    } else {
        "help".to_string()
    }
}

/// Determine the type of input
pub fn determine_type(input: &str) -> Result<InputType, Box<dyn std::error::Error>> {
    let file = File::open(input).unwrap();
    // Create a buffered reader for the file
    let reader = BufReader::new(file);

    // Read the first line of the file
    let mut lines = reader.lines();
    if let Some(line) = lines.next() {
        let first_line = line.unwrap();
        if first_line.starts_with('A') {
            return Ok(InputType::Segment);
        } else if first_line.starts_with('G') {
            // Count + and - in the firstline
            let mut pm = 0;
            for c in first_line.chars() {
                if c == '+' {
                    pm += 1;
                } else if c == '-' {
                    pm += 1;
                }
            }
            if pm == 2 {
                return Ok(InputType::Link);
            } else {
                return Ok(InputType::DirSegment);
            }
        } else if first_line.starts_with('S') {
            return Ok(InputType::Subgraph);
        } else {
            return Ok(InputType::Block);
        }
    } else {
        println!("The file is empty.");
    }
    Ok(InputType::Block)
}

pub fn read_file_lines(_file_path: &str, class: &InputType) -> std::io::Result<Vec<u64>> {
    let file = File::open("example.txt")?;

    // Create a buffered reader for efficient reading.
    let reader = BufReader::new(file);

    let mut vec_u64 = Vec::new();

    // Iterate over each line in the file.
    for line in reader.lines() {
        let line = line?[1..].to_string(); // Handle any errors in reading lines.
        let mut bools = [false; 2];
        let mut jj = [0; 2];
        let mut index = 0;
        for x in line.chars() {
            if x == '+' || x == '-' {
                if x == '+' {
                    bools[0] = true;
                }
                jj[0] = line[0..index].parse().unwrap();
                break;
            }
            index += 1;
        }
        if line.ends_with('+') || line.ends_with('-') {
            bools[1] = line.ends_with('+');
            jj[1] = line[index..line.len() - 2].parse().unwrap();
        } else {
            jj[1] = line[index..line.len() - 1].parse().unwrap();
        }

        if *class == InputType::Segment {
            vec_u64.push(jj[1] as u64);
        } else if *class == InputType::DirSegment {
            vec_u64.push(jj[0] as u64 * 2 + bools[0] as u64);
        } else if *class == InputType::Link {
            let u1 = jj[0] * 2 + bools[0] as u32;
            let u2 = jj[1] * 2 + bools[1] as u32;
            vec_u64.push(merge_u32_to_u64(u1, u2));
        }
    }

    Ok(vec_u64)
}

/// positional vector
///
/// question. where can I find a node in a path, which position
/// Slice of a vector
pub fn pos(graph: &Gfa<u32, (), ()>, name: String, class: &InputType) -> Vec<(u32, u64)> {
    for path in graph.paths.iter() {
        if name == path.name {
            let mut pos: u32 = 0;
            let mut vec_u64 = Vec::new();
            for i in 0..path.nodes.len() - 1 {
                pos += graph.get_sequence_by_id(&path.nodes[i]).len() as u32;
                let v1 = path.nodes[i];
                let v2 = path.dir[i];
                let v3 = path.nodes[i + 1];
                let v4 = path.dir[i + 1];

                if *class == InputType::Segment {
                    vec_u64.push((pos, v1 as u64));
                } else if *class == InputType::DirSegment {
                    vec_u64.push((pos, v1 as u64 * 2 + v2 as u64));
                } else if *class == InputType::Link {
                    let u1 = v1 * 2 + v2 as u32;
                    let u2 = v3 * 2 + v4 as u32;
                    vec_u64.push((pos, merge_u32_to_u64(u1, u2)));
                }
            }
            vec_u64.sort();
            return vec_u64;
        }
    }
    Vec::new()
}

/// Find the closest reference node
pub fn find_easy(
    graph: &Gfa<u32, (), ()>,
    class: &InputType,
    data_input: Vec<u64>,
    output: &str,
) -> Result<(), std::io::Error> {
    let file_out = File::create(output)?;
    let mut output_reader = std::io::BufWriter::new(file_out);

    for path in graph.paths.iter() {
        let pos1 = pos(graph, path.name.clone(), class);
        let mut i = 0;
        let mut j = 0;
        let mut last_node = 0;

        while i < data_input.len() && j < pos1.len() {
            if data_input[i] == pos1[j].1 {
                if pos1[j].0 != last_node {
                    last_node = pos1[j].0;
                    j += 1;
                    i += 1;
                } else {
                    j += 1;
                }
                writeln!(
                    output_reader,
                    "{}\t{}\t{}",
                    path.name, data_input[i], pos1[j].1
                )?;
            } else if data_input[i] < pos1[j].1 {
                i += 1;
            } else {
                j += 1;
            }
        }
    }
    Ok(())
}
