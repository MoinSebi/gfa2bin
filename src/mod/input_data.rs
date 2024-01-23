use crate::core::helper::{merge_u32_to_u64, Feature};

use std::fs::{File, read};
use std::io::{BufRead, BufReader};

pub struct FileData {
    pub data: Vec<u64>,
    pub feature: Feature,
}

impl FileData {
    pub fn new() -> Self {
        Self {
            data: Vec::new(),
            feature: Feature::Node,
        }
    }

    pub fn from_file(filename: &str) -> Self {
        let feature = get_type(filename);

        let mut data = Vec::new();
        let file = File::open(filename).expect("ERROR: CAN NOT READ FILE\n");
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line.unwrap();
            if feature == Feature::Edge {
                let ss = find_first_plus_minus(&line).unwrap();
                let s1 = &line[..ss];
                let s2 = &line[ss..ss + 1];
                let s3 = &line[ss + 1..line.len() - 1];
                let s4 = &line.chars().last().unwrap();
                let ss1 = s1.parse::<u64>().unwrap() * 2 + (s2 == "+") as u64;
                let ss2 = s3.parse::<u64>().unwrap() * 2 + (*s4 == '+') as u64;
                data.push(merge_u32_to_u64(ss1 as u32, ss2 as u32));
            } else if feature == Feature::DirNode {
                let s = line.ends_with('+');
                let s2 = line[..line.len() - 1].parse::<u64>().unwrap() * 2 + s as u64;
                data.push(s2)
            } else {
                data.push(line.parse::<u64>().unwrap());
            }
        }
        data.sort();
        Self { data, feature }
    }
}

fn find_first_plus_minus(input: &str) -> Option<usize> {
    input.chars().position(|c| c == '+' || c == '-')
}

pub fn get_type(file_path: &str) -> Feature {
    let file = File::open(file_path).expect("ERROR: CAN NOT READ FILE\n");

    // Parse plain text or gzipped file
    let reader = BufReader::new(file);

    // Read the first line of the file
    let first_line = reader.lines().next().unwrap().unwrap();
    let parts: Vec<&str> = first_line.split(|c| c == '+' || c == '-').collect();
    let last_letter = first_line.chars().last().unwrap();
    if last_letter == '+' || last_letter == '-' {
        if parts.len() == 1 {
            Feature::DirNode
        } else {
            Feature::Edge
        }
    } else {
        Feature::Node
    }
}

pub fn read_paths(filename: &str) -> Vec<String>{
    let mut paths = Vec::new();

    let file = File::open(filename).expect("ERROR: CAN NOT READ FILE\n");

    // Parse plain text or gzipped file
    let reader = BufReader::new(file);

    for x in reader.lines(){
        paths.push(x.unwrap());
    }
    return paths
}
