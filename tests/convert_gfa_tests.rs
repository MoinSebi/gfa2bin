use std::fs;
use assert_cmd::prelude::*; // Add methods on commands
use std::fs::File;
use std::io::Read;
use std::process::Command;

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -t e (type edge)
fn gfa_edges() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("data/example_data/testGraph2.gfa")
        .arg("-o")
        .arg("data/output/graph.edge")
        .arg("-f")
        .arg("edge")
        .arg("--pansn")
        .arg("#");
    cmd.assert().success();

    let mut b = File::open("data/output/graph.edge.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + (11 * 2));
    assert_eq!(buffer[3], 63);
    fs::remove_file("data/output/graph.edge.bed")?;
    fs::remove_file("data/output/graph.edge.bim")?;
    fs::remove_file("data/output/graph.edge.fam")?;


    Ok(())
}

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -f edge
fn gfa_nodes() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("/home/svorbrugg/code/gfa2bin/data/example_data/testGraph2.gfa")
        .arg("-o")
        .arg("data/output/graph.node")
        .arg("-f")
        .arg("node");
    cmd.assert().success();
    let mut b = File::open("data/output/graph.node.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + (8 * 2));
    assert_eq!(buffer[3], 127);
    assert_eq!(buffer[4], 0);
    fs::remove_file("data/output/graph.node.bed")?;
    fs::remove_file("data/output/graph.node.bim")?;
    fs::remove_file("data/output/graph.node.fam")?;

    Ok(())
}

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -f dirnode
fn gfa_dir_node() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("/home/svorbrugg/code/gfa2bin/data/example_data/testGraph2.gfa")
        .arg("-o")
        .arg("data/output/graph.dirnode")
        .arg("-f")
        .arg("dirnode");
    cmd.assert().success();
    let mut b = File::open("data/output/graph.dirnode.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + (8 * 2));
    assert_eq!(buffer[3], 127);
    assert_eq!(buffer[4], 0);

    fs::remove_file("data/output/graph.dirnode.bed")?;
    fs::remove_file("data/output/graph.dirnode.bim")?;
    fs::remove_file("data/output/graph.dirnode.fam")?;

    Ok(())
}

#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -f dirnode
fn gfa_dir_node_paths_remove() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("/home/svorbrugg/code/gfa2bin/data/example_data/testGraph2.gfa")
        .arg("-o")
        .arg("data/output/graph.dirnode.paths")
        .arg("-f")
        .arg("dirnode")
        .arg("--paths")
        .arg("/home/svorbrugg/code/gfa2bin/data/example_data/paths_paths.txt");
    cmd.assert().success();
    let mut b = File::open("data/output/graph.dirnode.paths.bed").unwrap();
    let mut buffer = Vec::new();
    b.read_to_end(&mut buffer).unwrap();
    assert_eq!(buffer.len(), 3 + 8);
    assert_eq!(buffer[3], 31);
    assert_eq!(buffer[4], 15);
    fs::remove_file("data/output/graph.dirnode.paths.bed")?;
    fs::remove_file("data/output/graph.dirnode.paths.bim")?;
    fs::remove_file("data/output/graph.dirnode.paths.fam")?;
    Ok(())
}



#[test]
/// Test for plink subcommand
/// -g (gfa)
/// -f dirnode
fn gfa_dir_node_bimbam() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("gfa2bin")?;
    cmd.arg("graph")
        .arg("-g")
        .arg("/home/svorbrugg/code/gfa2bin/data/example_data/testGraph2.gfa")
        .arg("-o")
        .arg("data/output/graph.dirnodes")
        .arg("-f")
        .arg("edge")
        .arg("--bimbam");
    cmd.assert().success();
    let mut b = File::open("data/output/graph.dirnodes.bimbam").unwrap();
    let mut content = String::new();
    b.read_to_string(&mut content).unwrap();
    assert_eq!(content.contains("1+2+"), true);
    let mut b = File::open("data/output/graph.dirnodes.pheno").unwrap();
    let mut content = String::new();
    b.read_to_string(&mut content).unwrap();
    assert_eq!(content.contains("a"), true);
    fs::remove_file("data/output/graph.dirnodes.bimbam")?;
    fs::remove_file("data/output/graph.dirnodes.pheno")?;

    Ok(())
}
