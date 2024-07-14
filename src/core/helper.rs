use crate::remove::input_data::find_first_plus_minus;
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
    Block,
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
            "block" => Feature::Block,
            _ => panic!("Not implemented"),
        }
    }

    pub fn string2u128(s: &str, feature: Feature, feature2: Option<Feature>) -> u128 {
        match feature {
            Feature::Node => s.parse().unwrap(),
            Feature::DirNode => {
                let s1 = s.ends_with('+');
                let s2 = s[..s.len() - 1].parse::<u128>().unwrap() * 2 + s1 as u128;
                s2
            }
            Feature::Edge => {
                let ss = find_first_plus_minus(s).unwrap();
                let s1 = &s[..ss];
                let s2 = &s[ss..ss + 1];
                let s3 = &s[ss + 1..s.len() - 1];
                let s4 = &s.chars().last().unwrap();
                let ss1 = s1.parse::<u128>().unwrap() * 2 + (s2 == "+") as u128;
                let ss2 = s3.parse::<u128>().unwrap() * 2 + (*s4 == '+') as u128;
                merge_u64_to_u128(ss1 as u64, ss2 as u64)
            }
            Feature::Alignment => 1,
            Feature::PWindow => 1,
            Feature::Block => 1,
            Feature::MWindow => 1,
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
            Feature::Block => "block".to_string(),
        }
    }

    /// Convert the "index"-u64 to a String
    pub fn to_string_u64(&self, input: u64) -> String {
        if *self == Feature::PWindow {
            let (left, right) = split_u64_to_u32s(input);

            "P".to_string() + &format_unsigned_as_string(left) + &format_unsigned_as_string(right)
        } else if *self == Feature::MWindow {
            let (left, right) = split_u64_to_u32s(input);
            return "M".to_string()
                + &format_unsigned_as_string(left)
                + &format_unsigned_as_string(right);
        } else if *self == Feature::Block {
            let (left, right) = split_u64_to_u32s(input);
            return "B".to_string()
                + &format_unsigned_as_string(left)
                + &format_unsigned_as_string(right);
        } else if Feature::Node == *self {
            input.to_owned().to_string()
        } else if *self == Feature::DirNode {
            return format_unsigned_as_string(input);
        } else if *self == Feature::Edge {
            let (left, right) = split_u64_to_u32s(input);

            return format_unsigned_as_string(left) + &format_unsigned_as_string(right);
        } else {
            let (left, right) = split_u64_to_u32s(input);
            return left.to_string() + "_" + right.to_string().as_str();
        }
    }

    pub fn identify_feature(pp: &str) -> (Feature, Option<Feature>) {
        let parts: Vec<&str> = pp
            .split(|c| c == '+' || c == '-')
            .filter(|s| !s.is_empty())
            .collect();
        let last_letter = pp.chars().last().unwrap();
        let parts2 = pp.split('_').collect::<Vec<&str>>();
        if pp.starts_with('P') {
            (Feature::PWindow, None)
        } else if pp.starts_with('M') {
            (
                Feature::MWindow,
                Some(Feature::identify_feature(pp[1..].split('_').next().unwrap()).0),
            )
        } else if pp.starts_with('B') {
            (Feature::Block, None)
        } else if parts2.len() == 2 {
            (Feature::Alignment, None)
        } else if last_letter == '+' || last_letter == '-' {
            if parts.len() == 1 {
                (Feature::DirNode, None)
            } else {
                (Feature::Edge, None)
            }
        } else {
            (Feature::Node, None)
        }
    }
}

pub fn read1(input: &str, f: Feature) -> (u64, u64) {
    match f {
        Feature::Node => (input.parse().unwrap(), 0),
        Feature::DirNode => {
            let last_char = &input[input.len() - 1..];
            let rest = &input[..input.len() - 1];

            (
                rest.parse::<u64>().unwrap() * 2 + (last_char == "+") as u64,
                0,
            )
        }
        Feature::Edge => {
            let ff = input.find(|c| c == '+' || c == '-').unwrap();
            let dir1 = &input[ff..ff];
            let dir2 = &input[input.len() - 1..input.len() - 1];

            let number1 = &input[..ff];
            let number2 = &input[ff + 1..];

            let numb1: u32 = number1.parse::<u32>().unwrap() * 2 + (dir1 == "+") as u32;
            let numb2: u32 = number2.parse::<u32>().unwrap() * 2 + (dir2 == "+") as u32;

            (merge_u32_to_u64(numb1, numb2), 0)
        }
        Feature::Alignment => {
            let ff: Vec<u32> = input.split('_').map(|a| a.parse().unwrap()).collect();

            (merge_u32_to_u64(ff[0], ff[1]), 0)
        }
        Feature::MWindow => {
            let ff: Vec<u32> = input.split('_').map(|a| a.parse().unwrap()).collect();
            (merge_u32_to_u64(ff[0], ff[2]), 0)
        }
        Feature::PWindow => {
            let ff: Vec<u32> = input.split('_').map(|a| a.parse().unwrap()).collect();
            (merge_u32_to_u64(ff[0], ff[2]), 0)
        }
        Feature::Block => {
            let ff: Vec<u32> = input.split('_').map(|a| a.parse().unwrap()).collect();
            (merge_u32_to_u64(ff[0], ff[2]), 0)
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

    result
}

pub fn merge_u32_to_u64(high: u32, low: u32) -> u64 {
    let high_u64 = u64::from(high);
    let low_u64 = u64::from(low);

    let result: u64 = (high_u64 << 32) | low_u64;

    result
}

pub fn merge_u64_to_u128(high: u64, low: u64) -> u128 {
    let high_u128 = (high as u128) << 64;
    let low_u128 = low as u128;
    high_u128 | low_u128
}

pub fn is_all_zeros(bitvector: &BitVec<u8, Lsb0>) -> bool {
    return bitvector.iter().all(|byte| !byte);
}

pub fn is_all_ones(bitvector: &BitVec<u8, Lsb0>) -> bool {
    return bitvector.iter().all(|byte| *byte);
}

pub fn average_vec_u16(vector: &[u16], rval: f32) -> f64 {
    let sum: u16 = vector.iter().sum();
    let len = vector.len() as f64;
    (sum as f64 / len) * rval as f64
}

pub fn percentile(data: &[u16], u: f64) -> f64 {
    let mut data2 = data.to_vec();
    data2.sort();

    let n = data.len() as f64;
    let p: usize = ((u / 100.0) * (n - 1.0)).floor() as usize;

    data2[p] as f64
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
        let median =
            (sorted_data[middle_index_1] as f64 + sorted_data[middle_index_2] as f64) / 2.0;
        median * rval as f64 / 100_f64
    } else {
        // If odd, return the middle element
        let middle_index = n / 2;
        sorted_data[middle_index] as f64 * rval as f64 / 100_f64
    }
}
