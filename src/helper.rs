pub fn binary2dec_bed(vecc: &[bool]) -> u8{
    let mut result: u8 = 0;
    let mut count = 0;
    for x in vecc.iter().rev(){
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

