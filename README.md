# gfa2bin

Convert a gfa to a plink file. This tool can also use a compressed coverage file from packing.  

## Installation
```
git clone https://github.com/MoinSebi/gfa2bin
cd gfa2bin
cargo build --release
./target/release/gfa2bin -h 
```


## Subcommands
### graph - Converting from graphs directly
Convert a graph in gfa format to a plink or bimbam format. You are able to specify which feature (```-f ```) you want to use as entry. We support nodes (1), edges (1+2+), and directed nodes (1+). Path can be merged to samples using the PanSN-spec, which is highly recommended.  
You are able to ignore certain path using the ```--path``` option. Nodes can also be ignored using the ```gfa2bin mod``` command after you have converted the graph into plink format.  

##### Workflow
We count the occurrence of each feature in each path/sample in the graph. This is done by iterating through the graph and counting the occurrence of each feature in each path. Dependent on output format, you are able to create a presence/absence matrix (plink) or a normalized matrix (bimbam). You can either provide an absolute value as threshold ```-a```, which will be used in all entries, or you provide a method (mean, median, percentile), which will be computed per row (feature). The resulting value will then be multiplied by the relative threshold. If no threshold is provided, you will receive a presence/absence plink or max-val normalized bimbam file. 
##### Diploid
We are able to provide information about ploidy. This is auto-detected when using the PanSN-spec. In PLINK files, ploidy can easily be represented by 11, 01, 10, 00. In a bimfile, we use the average of both "normalized" haplotypes.

##### BIMBAM
In bimbam format, the value used for normalization represents the 2.0 in the 0.0 to 2.0 range. As said above, we average the two haplotpyes after the "0-2" normalization. Not sure if this is right. 

##### Comment: 
In our experience there is no need for a columns (path) normalization, since samples/paths can contain extreme numbers of single features which mess up the normalization. If wanted, I can implement this in the feature.  

#### Example usage: 
````text
gfa2bin graph -g input.gfa -o output -f node --bimbam 
gfa2bin graph -g input.gfa -o output -f dirnode -m mean -r 50 --pansn '#'
````

### align - Using graph alignment

Convert coverage information from sequence to graph alignments to plink bed files. Either can use plain pack files directly (which will consume large amount of memory) or use one of the custom coverage file formats from packiong repository as input. The packing repository helps to reduce storage and can perform pre-processing on sample level.  
Comparable to graph subcommand, we offer additional normalization can be run when using value based input. THis normalization is then run on feature level (e.g. nodes or edges). 

#### Example usage: 


### Remove
Remove samples or entries from the plink files (bed, bim, fam). 

#### Example usage: 
````text
./target/release/gfa2bin remove -b input --samples samples.txt --genotypes genotypes_names.txt -o output_plink
````
Tip: Use gretl or any other tool to get a list of samples or entries with a specific statistic.

### Filter
Filter entries or samples from a plink file. 

**Samples can be filtered by:**
- Below certain sample number 
- Above certain sample number

**Genotypes can be filtered by:**
- MAF (major allele frequency)
- maf (minor allele frequency)

Genotypes can be filtered by:


### Split and merge
#### Split 
Split a single plink file (bed, bim, fam), into multiple parts of the same size. This might be prefered if the testing data set is very big and performing GWAS takes a lot of time and multiprocessing not is possible. 

#### Merge
Merging multiple plink files back together. Either from the above computation or any other splitting operation. Samples in all input files, must be in same order (similar fam order and names). Input is file of name of all bed files (fam and bim should have the same prefix).

### View 

Convert a plink bed file to a vcf-like file format. This method might be useful for general checking of the generated genotypes. File might be of huge size dependent on input. 

#### Example 

### Find 

Given a list of genotypes (e.g. significant nodes or edges) and graph, return the position (in bed format) of those paths, where such genotypes can be found. Each genotype will be listed as additional information in the bed file. If users might need more than just the exact position, additional --length information can be added, which will return in bigger intervals, adding the additinal lengthg to each site.  
The output is made for extracting the sequence from the initial sequence and blasting these back to a database to get more information about selected DNA segment (overlap with genes or other interesting regions). 

#### Example usage
````text
./target/release/gfa2bin nearest -g graph.gfa -p 'a' -o output.table.txt
````
#### Example output
| node | ref_node | distance | position | path     |
|------|----------|----------|----------|----------|
| 3    | 1        | 0        | 0        | a#1#Chr1 |
| 1    | 1        | -1       | 0        | a#1#Chr1 |
| 4    | 2        | 0        | 10       | a#1#Chr1 |
| 2    | 2        | -1       | 10       | a#1#Chr1 |
| 5    | 5        | -1       | 15       | a#1#Chr1 |
Comment: Distance is the distance between the input node and the reference node in base pairs. Position is the position of the reference node in the reference path. Distance of -1 means that the node is a reference node and 0 means that the node is one node away from the reference node (no nodes, bp in between).


### Nearest node 

Return the closest reference-node in resprect to the input node. A reference node is the clostest node which can be found on any given reference path. The result does additionally return reference position of this node, An example is shown shown below. 

#### Example usage 
````text
./target/release/gfa2bin nearest -g graph.gfa -p 'a' -o output.table.txt
````
#### Example output
| node | ref_node | distance | position | path     |
|------|----------|----------|----------|----------|
| 3    | 1        | 0        | 0        | a#1#Chr1 |
| 1    | 1        | -1       | 0        | a#1#Chr1 |
| 4    | 2        | 0        | 10       | a#1#Chr1 |
| 2    | 2        | -1       | 10       | a#1#Chr1 |
| 5    | 5        | -1       | 15       | a#1#Chr1 |
Comment: Distance is the distance between the input node and the reference node in base pairs. Position is the position of the reference node in the reference path. Distance of -1 means that the node is a reference node and 0 means that the node is one node away from the reference node (no nodes, bp in between). 


### Testing 
````text
cargo test
````