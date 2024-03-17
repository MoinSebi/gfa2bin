use bitvec::order::Lsb0;
use bitvec::prelude::BitVec;
use std::fmt::Display;

#[derive(Debug, Clone, Eq, PartialEq, Copy)]
pub enum Feature {
    Node,
    DirNode,
    Edge,
    Alignment,
    MWindow,
    PWindow,
}

impl Feature {
    pub fn from_str(s: &str) -> Self {
        match s {
            "node" => Feature::Node,
            "dirnode" => Feature::DirNode,
            "edge" => Feature::Edge,
            "alignment" => Feature::Alignment,
            "mwindow" => Feature::MWindow,
            "pwindow" => Feature::PWindow,
            _ => panic!("Not implemented"),
        }
    }

    pub fn to_string1(&self) -> String {
        match self {
            Feature::Node => "node".to_string(),
            Feature::DirNode => "dirnode".to_string(),
            Feature::Edge => "edge".to_string(),
            Feature::Alignment => "alignment".to_string(),
            Feature::MWindow => "mwindow".to_string(),
            Feature::PWindow => "pwindow".to_string(),
        }
    }
}

pub fn from_string(name_input: &str, ftype: Feature) -> u64 {
    if ftype == Feature::Node {
        name_input.parse().unwrap()
    } else if ftype == Feature::DirNode {
        let last_char = &name_input[name_input.len() - 1..];
        let rest = &name_input[..name_input.len() - 1];

        rest.parse::<u64>().unwrap() * 2 + (last_char == "+") as u64
    } else {
        let ff = name_input.find(|c| c == '+' || c == '-').unwrap();
        let dir1 = &name_input[ff..ff];
        let dir2 = &name_input[name_input.len() - 1..name_input.len() - 1];

        let number1 = &name_input[..ff];
        let number2 = &name_input[ff + 1..];

        let numb1: u32 = number1.parse::<u32>().unwrap() * 2 + (dir1 == "+") as u32;
        let numb2: u32 = number2.parse::<u32>().unwrap() * 2 + (dir2 == "+") as u32;

        merge_u32_to_u64(numb1, numb2)
    }
}

pub fn to_string1(input: u64, ftype: &Feature) -> String {
    if Feature::Node == *ftype {
        input.to_owned().to_string()
    } else if *ftype == Feature::DirNode {
        return format_unsigned_as_string(input);
    } else if *ftype == Feature::Edge {
        let (left, right) = split_u64_to_u32s(input);

        return format_unsigned_as_string(left) + &format_unsigned_as_string(right);
    } else {
        let (left, right) = split_u64_to_u32s(input);
        return left.to_string() + "_" + right.to_string().as_str();
    }
}

fn format_unsigned_as_string<T: Display + Into<u64>>(name: T) -> String {
    let name_u64 = name.into();
    format!(
        "{}{}",
        name_u64 / 2,
        if name_u64 % 2 == 1 { "+" } else { "-" }
    )
}

pub fn split_u64_to_u32s(value: u64) -> (u32, u32) {
    let low = value as u32;
    let high = (value >> 32) as u32;

    (high, low)
}

pub fn index2node_seq(vector: &[u32]) -> Vec<u64> {
    let mut result = Vec::new();
    let mut node_id = vector[0];
    let mut seq_count = 0;
    result.push(merge_u32_to_u64(node_id, seq_count));

    for &value in vector.iter().skip(1) {
        if value == node_id {
            seq_count += 1;
            result.push(merge_u32_to_u64(node_id, seq_count))
        } else {
            node_id = value;
            seq_count = 0;
            result.push(merge_u32_to_u64(node_id, seq_count))

        }
    }
    result.push(merge_u32_to_u64(node_id, seq_count));
    result
}

pub fn merge_u32_to_u64(high: u32, low: u32) -> u64 {
    let high_u64 = u64::from(high);
    let low_u64 = u64::from(low);

    let result: u64 = (high_u64 << 32) | low_u64;

    result
}

pub fn is_all_zeros(bitvector: &BitVec<u8, Lsb0>) -> bool {
    return bitvector.iter().all(|byte| !byte);
}

pub fn is_all_ones(bitvector: &BitVec<u8, Lsb0>) -> bool {
    return bitvector.iter().all(|byte| *byte);
}


pub fn wrapper_stats(vector: &[u16], method: &str, rval: u16) -> f64 {
    match method {
        "mean" => average_vec_u16(vector, rval),
        "median" => median(vector, rval),
        "percentile" => percentile(vector, rval as f64),
        _ => panic!("Method not implemented"),
    }
}



pub fn average_vec_u16(vector: &[u16], rval: u16) -> f64 {
    let sum: u16 = vector.iter().sum();
    let len = vector.len() as f64;
    (sum as f64 / len ) * rval as f64 / 100 as f64
}

pub fn percentile(data: &[u16], u: f64) -> f64 {

    let mut data2 = data.to_vec();
    data2.sort();

    let n = data.len() as f64;
    let p: usize = ((u / 100.0) * (n - 1.0)).floor() as usize;

    return data2[p] as f64;
}

pub fn median(data: &[u16], rval: u16) -> f64 {


    // Create a mutable copy of the data and sort it
    let mut sorted_data = data.to_vec();
    sorted_data.sort();

    let n = sorted_data.len();

    // Check if the number of elements is odd or even
    if n % 2 == 0 {
        // If even, return the average of the middle two elements
        let middle_index_1 = (n / 2) - 1;
        let middle_index_2 = n / 2;
        let median = (sorted_data[middle_index_1] as f64 + sorted_data[middle_index_2] as f64) / 2.0;
        median * rval as f64 / 100 as f64
    } else {
        // If odd, return the middle element
        let middle_index = n / 2;
        sorted_data[middle_index] as f64 * rval as f64 / 100 as f64
    }
}



