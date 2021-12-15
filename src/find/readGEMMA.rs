use std::fs::File;
use std::io::{BufReader, BufRead};

fn readTable(filename: &str, col_pval: usize, col_pos: usize) -> (Vec<f64>, Vec<usize>){
    let file = File::open(filename).expect("ERROR: CAN NOT READ FILE\n");
    let reader = BufReader::new(file);
    let mut vector_pval:Vec<f64> = Vec::new();
    let mut vector_pos: Vec<usize> = Vec::new();
    for (i, line) in reader.lines().enumerate(){
        let l = line.unwrap();
        if i != 0{
            let line_split: Vec<&str> = l.split("\t").collect();
            let pval: f64 = line_split[col_pval].parse().unwrap();
            vector_pval.push(pval);

            let pos: f64 = line_split[col_pos].parse().unwrap();
            vector_pval.push(pos);
        }
    }
    return (vector_pval, vector_pos)
}