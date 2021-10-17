
use std::fs::File;
use std::io::Read;
use std::fs;
use std::convert::TryInto;

pub struct read_in2{
    pub ty: String,
    pub name: String,
    pub cc: Vec<bool>,
}

pub struct read_in3{
    pub ty: String,
    pub name: String,
    pub cc: Vec<u16>,
}

/// get shit byte by byte
///
pub fn get_file_as_byte_vec(filename: &str) -> Vec<u8> {
    let mut f = File::open(&filename).expect("no file found");
    let metadata = fs::metadata(&filename).expect("unable to read metadata");
    let mut buffer = vec![0; metadata.len() as usize];
    f.read(&mut buffer).expect("buffer overflow");


    buffer
}

pub fn wrapper2(buffer: &Vec<u8>) -> Vec<read_in2>{
    // total length 73 + len
    let length = read_be_u32(&mut &buffer[3..7]);
    println!("{}", length);
    let oo = buffer.chunks((length + 73) as usize );
    println!("How many samples: {}", oo.len());
    let mut jo: Vec<read_in2> = Vec::new();
    for x in oo.into_iter(){
        let u = get_meta(x);
        let c = get_bin(x);
        jo.push(read_in2{name: u.3, ty: "dunno".to_string(), cc: c});
    }
    return jo

}

pub fn wrapper3(buffer: &Vec<u8>) -> Vec<read_in3>{
    // total length 73 + len
    let length = read_be_u32(&mut &buffer[3..7]);
    println!("{}", length);
    let oo = buffer.chunks((length + 73) as usize );
    println!("How many samples: {}", oo.len());
    let mut jo: Vec<read_in3> = Vec::new();
    for x in oo.into_iter(){
        let u = get_meta(x);
        let c = get_u16(x);
        jo.push(read_in3{name: u.3, ty: "dunno".to_string(), cc: c});
    }
    return jo

}

/// get the meta data
/// see definition
pub fn get_meta(buffer: & [u8]) -> (bool, u32, u16, String){
    let cov = buffer[3];
    let length = read_be_u32(&mut &buffer[3..7]);
    let thresh = read_be_u16(&mut &buffer[7..9]);
    let name = byte_to_string(&mut &buffer[9..73]);

    (cov == 1, length, thresh, name)

}

pub fn get_bin(buffer: & [u8]) -> Vec<bool>{
    let mut j: Vec<bool> = Vec::new();
    for x in buffer[73..].iter(){
        j.extend(byte_to_bitvec(&x));
    }
    j
}

pub fn get_u16(buffer: & [u8]) -> Vec<u16>{
    let mut j: Vec<u16> = Vec::new();
    let g = buffer[73..].chunks(2);

    for x in g.into_iter(){

        j.push(test1(& x));
    }
    j
}

fn byte_to_bitvec(buf: &u8) -> Vec<bool>{
    let mut h: Vec<bool> = Vec::new();
    let mut n = buf.clone();
    while (n > 0){
        h.push((n%2)  == 1);
        n = n/2
    }
    for x in 0..(8-h.len()){
        h.insert(0, false);
    }
    h
}


fn byte_to_string(input: &[u8]) -> String {
    let mut o = "".to_string();
    for x in input.iter(){

        o.push(x.clone() as char);

    }
    return o
}




fn read_be_u32(input: & mut &[u8]) -> u32 {
    let (int_bytes, rest) = input.split_at(std::mem::size_of::<u32>());
    *input = rest;
    u32::from_be_bytes(int_bytes.try_into().unwrap())
}

pub fn read_be_u16(input: &mut &[u8]) -> u16 {
    let (int_bytes, rest) = input.split_at(std::mem::size_of::<u16>());
    *input = rest;
    u16::from_be_bytes(int_bytes.try_into().unwrap())
}

fn test1(vector: &[u8]) -> u16{
    let number = ((vector[0] as u16) << 8) | vector[1] as u16;
    number
}


#[cfg(test)]
mod tests {
    use crate::pack_reader::{get_file_as_byte_vec, wrapper2};
    use crate::matrix_wrapper::matrix_node_coverage;

    #[test]
    fn it_works() {
        let h = get_file_as_byte_vec("/home/svorbrugg_local/Rust/packing/pack.2.bin");
        let h = get_file_as_byte_vec("/home/svorbrugg_local/Rust/packing/test.cov2.bin");
        let h2 = wrapper2(&h);
        let mut h3 = matrix_node_coverage(h2);
        h3.write("bed", "test.1.2", "node");

    }
}
