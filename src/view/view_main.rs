use std::fs::File;
use clap::ArgMatches;
use crate::core::core::MatrixWrapper;
use std::io::Write;
use bitvec::vec::BitVec;
use log::info;

/// View main function
///
/// Convert a bed file to VCF file
pub fn view_main(matches: &ArgMatches) {
    let plink_file = matches.value_of("plink").unwrap();
    let output_prefix = matches.value_of("output").unwrap();
    // Read the bed file


    info!("Reading Plink file");
    let mut mw = MatrixWrapper::new();
    let feature = mw.feature;
    mw.bfile_wrapper(plink_file);

    info!("Writing output (vcf)");
    mw.write_vcf(output_prefix);



}

impl MatrixWrapper{


    /// Write a VCF file with dumb header
    pub fn write_vcf(&self, filename: &str) {
        let file = File::create(filename).unwrap();
        let mut writer = std::io::BufWriter::new(file);
        write!(writer, "##fileformat=VCFv4.2\n");
        write!(writer, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");

        let fam_entries = self.fam_entries.iter().map(|x| x.split_whitespace().next().unwrap().to_string()).collect::<Vec<String>>();
        let bims = self.bim_entries.iter().map(|x| x.split_whitespace().nth(3).unwrap().to_string()).collect::<Vec<String>>();

        write!(writer, "{}", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t".to_string() + &fam_entries.join("\t") + "\n");

        for (x, y) in self.matrix_bit.iter().zip(bims.iter()){
            write!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", "graph", y, ".", y,"-",  self.feature.to_string1(), "PASS", "GT=".to_string() + &(x.len()/2).to_string(), "GT", bitvec2vcf_string(x)).expect("Error writing to file");
        }
    }
}


// Convert bitvec to vcf string
pub fn bitvec2vcf_string(bitvec: &BitVec<u8>) -> String {
    let mut s = String::new();
    let chunks = bitvec.chunks(2);
    for x in chunks {
        if x[0] == false && x[1] == false  {
            s.push_str("0/0");
        } else if x[0] == true && x[1] == false {
            s.push_str("0/1");
        } else if x[0] == false && x[1] == true  {
            s.push_str("1/0");
        } else if x[0] == true && x[1] == true {
            s.push_str("1/1");
        }
        s.push_str("\t");
    }
    s
}
