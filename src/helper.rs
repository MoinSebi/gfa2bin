use std::fs::File;
use std::io::{Write, BufWriter};
use bitvec::bitvec;
use bitvec::macros::internal::funty::Fundamental;
use bitvec::order::Lsb0;
use bitvec::vec::BitVec;
use log::info;
use packing_lib::helper::u8_u16;
use packing_lib::reader::get_file_as_byte_vec;


pub fn binary2dec_bed(vecc: &[bool]) -> u8{
    let mut result: u8 = 0;
    let mut count = 0;
    for x in vecc.iter(){
        let t: u8 = 2;
        result += (t.pow(count as u32)) * (*x as u8);
        count += 1;
        result += (t.pow(count as u32)) * (*x as u8);
        count += 1;
    }
    result
}

#[allow(dead_code)]
pub fn tuple2string(t: (u32, bool)) -> String{
    let j: String =  t.0.to_string() + &bool2string(t.1);
    j
}

#[allow(dead_code)]
pub fn bool2string(b: bool) -> String{
    if b{
        return "+".to_string();
    }
    else {
        return "-".to_string();
    }
}


#[allow(dead_code)]
pub fn transpose<T>(v: &Vec<Vec<T>>) -> Vec<Vec<T>>
    where
        T: Clone,
{
    assert!(!v.is_empty());
    (0..v[0].len())
        .map(|i| v.iter().map(|inner| inner[i].clone()).collect::<Vec<T>>())
        .collect()
}

pub fn trans4(v: &Vec<bitvec::vec::BitVec>){
    let mut h: Vec<bitvec::vec::BitVec> = Vec::new();
    for x in v.iter(){
        let mut o10: BitVec<u8, Lsb0> = BitVec::new();
        for y in x.iter(){
            o10.push(y.as_bool());
        }
    }

}
pub fn trans2<T>(v: &Vec<Vec<T>>) -> Vec<Vec<T>>
where
T: Clone
{
    info!("Transposing");
    let mut o: Vec<Vec<T>> = Vec::new();
    for x in 0..v[0].len(){
        let mut o2: Vec<T> = Vec::new();
        for y in 0..v.len(){
            o2.push(v[y][x].clone());
        }
        o.push(o2);
    }
    return o;
}

pub fn trans3<T>(v: &[Vec<T>]) -> Vec<Vec<T>>
    where
        T: Clone
{
    info!("Transposing");
    let mut o: Vec<Vec<T>> = Vec::new();
    for x in 0..v[0].len(){
        let mut o2: Vec<T> = Vec::new();
        for y in 0..v.len(){
            o2.push(v[y][x].clone());
        }
        o.push(o2);
    }
    return o;
}


pub fn trans5(v: &Vec<bitvec::vec::BitVec>) -> Vec<bitvec::vec::BitVec>{
    info!("Transposing");
    let mut o: Vec<bitvec::vec::BitVec> = Vec::new();

    for x in 0..v[0].len(){
        let mut o2: bitvec::vec::BitVec = bitvec::vec::BitVec::new();
        for y in 0..v.len(){
            o2.push(v[y][x].clone());
        }
        o.push(o2);
    }
    return o;
}



pub fn get_thresh(filename: &str) -> u16{
    let size = u8_u16(&mut & get_file_as_byte_vec(filename)[7..9]);
    size
}

