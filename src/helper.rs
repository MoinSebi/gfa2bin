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
