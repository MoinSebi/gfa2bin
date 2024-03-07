use crate::core::core::MatrixWrapper;
use crate::core::helper::{from_string, Feature};
use bitvec::order::Lsb0;
use bitvec::prelude::BitVec;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

impl MatrixWrapper {
    pub fn bfile_wrapper(&mut self, filename: &str) {
        self.read_fam(&format!("{}{}", filename, ".fam"));
        self.read_bim(&format!("{}{}", filename, ".bim"));
        self.read_bed(&format!("{}{}", filename, ".bed"));
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
    pub fn read_bed(&mut self, filename: &str) {
        let mut file = File::open(filename).expect("no file found");
        let metadata = fs::metadata(filename).expect("unable to read metadata");
        let mut buffer = vec![0; metadata.len() as usize];

        file.read_exact(&mut buffer).expect("buffer overflow");
        // cut off first 3 bytes
        buffer = buffer[3..].to_vec();

        // do i need this
        let mut num = self.sample_names.len() / 4;
        if (self.sample_names.len() % 4) > 0 {
            num += 1;
        }
        // Each chunk
        let chunks = buffer.chunks(num);
        for chunk in chunks.into_iter() {
            let mut bv: BitVec<u8, Lsb0> = BitVec::from_slice(chunk);
            bv.truncate(self.sample_names.len() * 2);

            self.matrix_bit.push(bv);
        }
    }

    /// Read plink fam file
    ///
    /// More information about the file format: https://www.cog-genomics.org/plink/1.9/formats#fam
    /// In general [FamilyID, withinFamId (wfi), wfi or father, wfi or mother, sex_code, phenotype value]
    /// Example: tair10 0   0   0   0   0 321.32
    /// 0 = unknown or not in data set
    ///
    ///
    pub fn read_fam(&mut self, filename: &str) {
        let buf = BufReader::new(File::open(filename).expect("Unable to open file"));
        // Split the file at the newlines
        for line in buf.lines() {
            let vv = line.unwrap();
            let splits = vv.split_whitespace().collect::<Vec<&str>>();
            self.fam_entries.push(splits.join("\t"));
            self.sample_names.push(splits[0].to_string());
        }
    }

    /// Read plink BIM file
    ///
    /// More information about the file: https://www.cog-genomics.org/plink/1.9/formats#bim
    /// Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
    /// Variant identifier
    /// Position in morgans or centimorgans (safe to use dummy value of '0')
    /// Base-pair coordinate (1-based; limited to 231-2)
    /// Allele 1 (corresponding to clear bits in .bed; usually minor)
    /// Allele 2 (corresponding to set bits in .bed; usually major)
    /// [Chromosome code, variant identifier, position in morgans, bp, allele1, allele2]
    pub fn read_bim(&mut self, filename: &str) {
        let ff = get_type_bim(filename);

        let file = BufReader::new(File::open(filename).expect("Unable to open file"));

        for (_i, x) in file.lines().enumerate() {
            let data = x.unwrap();
            let a = from_string(data.split_whitespace().nth(3).unwrap(), ff);
            self.geno_names.push(a);
            self.bim_entries.push(data);
        }
    }
}

pub fn get_type_bim(file_path: &str) -> Feature {
    let file = File::open(file_path).expect("ERROR: CAN NOT READ FILE\n");

    // Parse plain text or gzipped file
    let reader = BufReader::new(file);

    // Read the first line of the file
    let first_line = reader.lines().next().unwrap().unwrap();
    let pp = first_line.split_whitespace().nth(3).unwrap();
    let parts: Vec<&str> = pp.split(|c| c == '+' || c == '-').collect();
    let last_letter = pp.chars().last().unwrap();
    if last_letter == '+' || last_letter == '-' {
        if parts.len() == 1 {
            Feature::DirNode
        } else {
            Feature::Edge
        }
    } else {
        Feature::Node
    }
}
