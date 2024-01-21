use byteorder::{BigEndian, ByteOrder};

use packing_lib::reader::get_file_as_byte_vec;

#[allow(dead_code)]
pub fn binary2dec_bed(vecc: &[bool]) -> u8 {
    let mut result: u8 = 0;
    let mut count = 0;
    for x in vecc.iter() {
        let t: u8 = 2;
        result += (t.pow(count as u32)) * (*x as u8);
        count += 1;
        result += (t.pow(count as u32)) * (*x as u8);
        count += 1;
    }
    result
}

pub fn make_dir_name(maxval: &usize) -> Vec<(usize, bool)> {
    let mut f = Vec::new();
    for x in 0..*maxval {
        f.push((x, false));
        f.push((x, true));
    }
    f
}

pub fn get_thresh(filename: &str) -> u16 {
    BigEndian::read_u16(&get_file_as_byte_vec(filename)[7..9])
}
