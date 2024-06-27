use crate::core::helper::Feature;

use std::fs::File;
use std::io::{BufRead, BufReader};

pub struct FileData {
    pub data: Vec<u64>,
    pub window: Vec<u32>,
    pub feature: Feature,
}

impl FileData {
    pub fn from_file(filename: &str) -> Self {
        let feature = get_type(filename);

        let mut data = Vec::new();
        let mut data3 = Vec::new();
        let mut data2 = Vec::new();
        let file = File::open(filename).expect("ERROR: CAN NOT READ FILE\n");
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let a = Feature::string2u64(&line.unwrap(), feature.0, feature.1);
            data.push((a.0, a.1));
        }

        // Sort very important
        data.sort();

        if feature.0 == Feature::MWindow || feature.0 == Feature::Block {
            for x in data.iter() {
                data2.push(x.1);
                data3.push(x.0);
            }
        } else {
            data3 = data.into_iter().map(|x| x.0).collect();
        }
        Self {
            data: data3,
            window: data2,
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
