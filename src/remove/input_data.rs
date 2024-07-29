use crate::core::helper::Feature;

use std::fs::File;
use std::io::{BufRead, BufReader};

pub struct FileData {
    pub data: Vec<u128>,
    pub feature: Feature,
}

impl FileData {
    pub fn from_file(filename: &str) -> Self {
        let feature = get_type(filename);

        let mut data = Vec::new();
        let data3 = Vec::new();
        let file = File::open(filename).expect("ERROR: CAN NOT READ FILE\n");
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let a = Feature::string2u128(&line.unwrap(), feature.0, feature.1);
            data.push(a);
        }

        // Sort very important
        data.sort();

        Self {
            data: data3,
            feature: feature.0,
        }
    }
}

pub fn find_first_plus_minus(input: &str) -> Option<usize> {
    input.chars().position(|c| c == '+' || c == '-')
}

pub fn get_type(file_path: &str) -> (Feature, Option<Feature>) {
    let file = File::open(file_path).expect("ERROR: CAN NOT READ FILE\n");

    // Parse plain text or gzipped file
    let reader = BufReader::new(file);

    // Read the first line of the file
    let first_line = reader.lines().next().unwrap().unwrap();
    Feature::identify_feature(&first_line)
}

pub fn read_paths(filename: &str) -> Vec<String> {
    let mut paths = Vec::new();

    let file = File::open(filename).expect("ERROR: CAN NOT READ FILE\n");

    // Parse plain text or gzipped file
    let reader = BufReader::new(file);

    for x in reader.lines() {
        paths.push(x.unwrap());
    }
    paths
}
