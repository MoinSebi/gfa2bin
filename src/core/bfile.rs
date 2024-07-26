use crate::core::core::MatrixWrapper;
use crate::core::helper::Feature;
use bitvec::order::Lsb0;
use bitvec::prelude::BitVec;
use gfa_reader::Pansn;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read};
use std::{fs, io};

/// Read number of lines
pub fn count_lines(file_path: &str) -> Result<usize, std::io::Error> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    // Count the lines using iterator folding
    let num_lines = reader.lines().count();

    Ok(num_lines)
}

impl MatrixWrapper {
    pub fn bfile_wrapper(&mut self, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut bim_count = count_lines(&format!("{}{}", filename, ".bim"))?;
        let mut fam_count = count_lines(&format!("{}{}", filename, ".fam"))?;

        self.read_bed(&format!("{}{}", filename, ".bed"), 1, 1)?;
        Ok(())
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
    pub fn read_bed(
        &mut self,
        filename: &str,
        samples_number: usize,
        snp_number: usize,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::open(filename)?;
        let metadata = fs::metadata(filename)?;
        let mut buffer = vec![0; metadata.len() as usize];

        file.read_exact(&mut buffer)?;
        // cut off first 3 bytes
        buffer = buffer[3..].to_vec();

        // do i need this
        let mut num = samples_number / 4;
        if (samples_number % 4) > 0 {
            num += 1;
        }
        // Each chunk
        let chunks = buffer.chunks(num);
        for chunk in chunks.into_iter() {
            let mut bv: BitVec<u8, Lsb0> = BitVec::from_slice(chunk);
            bv.truncate(samples_number * 2);

            self.matrix_bit.push(bv);
        }
        assert_eq!(self.matrix_bit.len(), snp_number);
        Ok(())
    }

    /// Read number of lines
    fn count_lines(file_path: &str) -> Result<usize, std::io::Error> {
        let file = File::open(file_path)?;
        let reader = BufReader::new(file);

        // Count the lines using iterator folding
        let num_lines = reader.lines().count();

        Ok(num_lines)
    }
}

pub fn get_type_bim(file_path: &str) -> (Feature, Option<Feature>) {
    let file = File::open(file_path).expect("ERROR: CAN NOT READ FILE\n");

    // Parse plain text or gzipped file
    let reader = BufReader::new(file);

    // Read the first line of the file
    let first_line = reader.lines().next().unwrap().unwrap();
    let first_line = first_line.split_whitespace().nth(3).unwrap();
    Feature::identify_feature(first_line)
}

use std::io::Write;
pub fn write_dummy_fam(pansn: &Pansn<u32, (), ()>, outfile: &str) -> Result<(), io::Error> {
    let mut file = File::create(outfile)?;
    let mut bufwriter = BufWriter::new(file);
    for (i, x) in pansn.genomes.iter().enumerate() {
        writeln!(bufwriter, "{}\t{}\t0\t0\t0\t-9", i, x.name)?;
    }
    Ok(())
}
