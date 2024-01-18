use std::fs;
use std::fs::File;
use std::io::{BufReader, BufRead, Read};
use bitvec::order::Msb0;
use bitvec::prelude::BitVec;
use crate::core::helper::GenoName;
use crate::MatrixWrapper;


/// Wrapper for reading all entries of a plink genotype file (bed, bam, fam)
/// matrix_bin -> bim, bam
/// column_name
/// helper
pub fn bfile_wrapper(filename: &str, matrix: &mut MatrixWrapper, names: &mut Vec<String>) {
    read_fam(&format!("{}{}", filename, ".fam"), matrix);
    read_bim(&format!("{}{}", filename, ".bim"), names);
    read_bed(&format!("{}{}", filename, ".bed"), matrix,  names.len());
}

/// Read plink fam file
///
/// More information about the file format: https://www.cog-genomics.org/plink/1.9/formats#fam
/// In general [FamilyID, withinFamId (wfi), wfi or father, wfi or mother, sex_code, phenotype value]
/// Example: tair10 0   0   0   0   0 321.32
/// 0 = unknown or not in data set
///
///
pub fn read_fam(filename: &str, matrix: &mut MatrixWrapper){
    let data = fs::read_to_string(filename).expect("Unable to read file");

    // Split the file at the newlines
    let lines_vec: Vec<_> = data.split("\n").collect();
    println!("GFA2BIN: Fam file has {} entries.", lines_vec.len());

    // Convert to String Vector
    let lines_vec_string: Vec<_> = lines_vec.into_iter().map(|e| e.to_string()).collect();

    // Add name to thing
    for (i, x) in lines_vec_string.into_iter().enumerate(){
        let line_split: Vec<_> = x.split("\t").collect();
        matrix.geno_map.insert(GenoName{name: i as u64}, i);
    }
}

/// Read plink BIM file
///
/// More information about the file: https://www.cog-genomics.org/plink/1.9/formats#bim
/// [Chromosome code, variant identifier, position in morgans, bp, allele1, allele2]
/// allele code can hold more than one identifier
pub fn read_bim(filename: &str, snp_names: &mut Vec<String>){
    let file = BufReader::new(File::open(filename).expect("Unable to open file"));


    for x in file.lines() {
        let f: Vec<String> = x.unwrap().split("\t").map(| s | s.to_string()).collect();
        let ff = f[0].clone();
        let ff2 = f[3].clone();

        snp_names.push(format!("{}_{}", ff.to_string(), ff2));
    }
    println!("dasd{}", snp_names.len());
}



/// Read plink bed file
///
/// More information about the file format: https://www.cog-genomics.org/plink/1.9/formats#bed
/// Additional information: https://zzz.bwh.harvard.edu/plink/binary.shtml
///
/// File starts with:  01101100 00011011 00000001
/// Then genotype data: 00 homo (ref), 01 hetero, 11 homo (alternative), 10 missing
///
///      01101100
///      HGFEDCBA
///
///            AB   00  -- homozygote (first)
///          CD     11  -- other homozygote (second)
///        EF       01  -- heterozygote (third)
///      GH         10  -- missing genotype (fourth)
pub fn read_bed(filename: &str, matrix_w: & mut MatrixWrapper, numbsnp: usize) {

    let mut file = File::open(filename).expect("no file found");
    let metadata = fs::metadata(&filename).expect("unable to read metadata");
    let mut buffer = vec![0; metadata.len() as usize];

    file.read_exact(&mut buffer).expect("buffer overflow");
    // cut off first 3 bytes
    buffer = buffer[3..].to_vec();

    // do i need this
    let mut num = matrix_w.geno_map.len()/4;
    if (matrix_w.geno_map.len() % 4) > 0{
        num += 1;
    }
    // Each chunk
    let chunks = buffer.chunks(num);

    println!("{}", numbsnp);
    for chunk in chunks.into_iter() {
        let bv: BitVec<u8, Msb0> = BitVec::from_slice(&chunk[..]);
        let mut dd: BitVec<u8, Msb0> = BitVec::new();

        for (i, x) in bv.iter().step_by(2).enumerate(){
            if i < num as usize {
                dd.push(*x);
            }

        }
        matrix_w.matrix_bin.push(dd);
    }
}

