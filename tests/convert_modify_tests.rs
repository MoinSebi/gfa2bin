use assert_cmd::Command;
use predicates::prelude::predicate;

#[test]
/// Check if the normal convert command it working
fn pack_filter() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("convert")
        .arg("-p")
        .arg("/home/svorbrugg_local/Rust/packing/testing/jo3.bin.zst")
        .arg("-o")
        .arg("data/output/pack2")
        .arg("--filter");
    cmd.assert().success();
    //fs::remove_file("example_data/test3.bubble.stats")?;
    //fs::remove_file("example_data/test3.bubble.txt")?;
    //fs::remove_file("example_data/test3.traversal.bed")?;

    Ok(())
}