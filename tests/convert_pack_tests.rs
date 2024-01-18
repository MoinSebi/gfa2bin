use assert_cmd::Command;
use predicates::prelude::predicate;

#[test]
/// Check if the normal plink command it working
fn pack_unknown() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("plink")
        .arg("-p")
        .arg("/home/svorbrugg_local/Rust/packing/testing/jo1.bin.zst")
        .arg("-o")
        .arg("data/output/pack1");
    cmd.assert().success();
    //cmd.assert().stderr(predicate::str::contains("No file with such name"));
    //fs::remove_file("example_data/test3.bubble.stats")?;
    //fs::remove_file("example_data/test3.bubble.txt")?;
    //fs::remove_file("example_data/test3.traversal.bed")?;

    Ok(())
}

#[test]
/// Check if the normal plink command it working
fn pack_unknown2() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("plink")
        .arg("-p")
        .arg("/home/svorbrugg_local/Rust/packing/testing/jo3.bin.zst")
        .arg("-o")
        .arg("data/output/pack2")
        .arg("--split")
        .arg("5");
    cmd.assert().success();
    //cmd.assert().stderr(predicate::str::contains("No file with such name"));
    //fs::remove_file("example_data/test3.bubble.stats")?;
    //fs::remove_file("example_data/test3.bubble.txt")?;
    //fs::remove_file("example_data/test3.traversal.bed")?;

    Ok(())
}
