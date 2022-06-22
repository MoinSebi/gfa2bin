use std::fs::File;
use std::io::{BufRead, BufReader};

pub struct  Fam{
    pub family_id: Vec<String>,
}

impl Fam {
    pub fn from_file(filename: &str) -> Self{
        let file = File::open(filename).expect("ERROR: CAN NOT READ FILE\n");
        let reader = BufReader::new(file);
        let mut lines = vec![];
        for line in reader.lines(){
            lines.push(line.unwrap());
        }
        let k: Vec<String> = lines.iter().map(|x| x.split("\t").nth(0).unwrap().to_string()).collect();
        Self{
            family_id: k,
        }
    }

    pub fn new() -> Self{
        Self{
            family_id: Vec::new(),
        }
    }
}