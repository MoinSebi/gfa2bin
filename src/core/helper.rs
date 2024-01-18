use std::fmt::Display;

#[derive(Debug, Clone, Eq, PartialEq, Copy)]
pub enum Feature {
    Node,
    DirNode,
    Edge,

}

impl Feature {
    pub fn from_str(s: &str) -> Self {
        match s {
            "node" => Feature::Node,
            "dirnode" => Feature::DirNode,
            "edge" => Feature::Edge,
            _ => panic!("Not implemented"),
        }
    }

    pub fn to_string(&self) -> String{
        match self {
            Feature::Node => "node".to_string(),
            Feature::DirNode => "dirnode".to_string(),
            Feature::Edge => "edge".to_string(),
        }
    }
}

#[derive(Eq, Hash, PartialEq, Debug, Clone)]
pub struct GenoName {
    pub name: u64
}

impl GenoName {

    pub fn new() -> Self{
        return GenoName {name: 0};

    }
    pub fn from_string(name_input: &str, ftype: Feature) -> Self{
        if ftype == Feature::Node{
            GenoName {
                name: name_input.parse().unwrap()
            }
        } else if ftype == Feature::DirNode{
            let last_char = &name_input[name_input.len()-1..];
            let rest = &name_input[..name_input.len()-1];


            GenoName {
                name: rest.parse::<u64>().unwrap()*2 + (last_char == "+") as u64
            }
        } else {
            let ff = name_input.find(|c| c == '+' || c == '-').unwrap();
            let dir1 = &name_input[ff..ff];
            let dir2 = &name_input[name_input.len()-1..name_input.len()-1];

            let number1 = &name_input[..ff];
            let number2 = &name_input[ff+1..];

            let numb1: u32 = number1.parse::<u32>().unwrap()*2 + (dir1 == "+") as u32;
            let numb2: u32 = number2.parse::<u32>().unwrap()*2 + (dir2 == "+") as u32;
            GenoName {
                name: merge_u32_to_u64(numb1, numb2)
            }
        }
    }


    pub fn to_string(&self, ftype: &Feature) -> String{
        if Feature::Node == *ftype{
            return self.name.to_string()
        } else if *ftype == Feature::DirNode{
            return format_unsigned_as_string(self.name)
        } else {
            let (left, right) = split_u64_to_u32s(self.name);

            return format_unsigned_as_string(left) + &format_unsigned_as_string(right)
        }
    }
}

fn format_unsigned_as_string<T: Display + Into<u64>>(name: T) -> String {
    let name_u64 = name.into();
    format!("{}{}", name_u64 / 2, if name_u64 % 2 == 1 { "-" } else { "+" })
}


fn split_u64_to_u32s(value: u64) -> (u32, u32) {
    let low = value as u32;
    let high = (value >> 32) as u32;

    (high, low)
}

fn merge_u32_to_u64(high: u32, low: u32) -> u64 {
    let high_u64 = u64::from(high);
    let low_u64 = u64::from(low);

    let result: u64 = (high_u64 << 32) | low_u64;

    result
}