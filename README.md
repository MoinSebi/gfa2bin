# gfa2bin

Convert a gfa to a plink file. This tool can also use a compressed coverage file from packing.  


## Subcommands
### graph - Converting from graphs directly
Convert a graph in gfa format to a plink or bimbam format. You are able to specify which feature (```-f ```) you want to use as feature. We support nodes (1), edges (1+2+), and directed nodes (1+). Path can be merged to samples using the PanSN-spec, which is highly recommended.  
You are able to ignore certain path using the ```--path``` option. Nodes can also be ignored using the ```gfa2bin mod``` command after you have converted the graph into plink format.  

##### Workflow
We count the occurence of each feature in each path/sample in the graph. This is done by iterating through the graph and counting the occurence of each feature in each path. Dependent on output format, you are able to create a presence/absence matrix (plink) or a a normalized matrix (bimbam). You can either provide a absolute value as threshold, which will be used in all entries, or you provide a method (mean, median, percentile, maximum), which will be computed feature-wise. The resulting value will then be multiplied by the relative threshold. If no threshold is provided, you will receive a presence/absence plink or max-val normalized bimbam file. 
##### Diploid
We are able to provide information about ploidy. This is auto-detected when using the PanSN-spec. In PLINK files, ploidy can easily be represented by 11, 01, 10, 00. In a bimfile, we use the average of both "normalized" haplotypes.

##### BIMBAM
In bimbam format, the value used for normalization represents the 2.0 in the 0.0 to 2.0 range. As said above, we average the two haplotpyes after the "0-2" normalization. Not sure if this is right. 



#### Example usage: 
````text
gfa2bin graph -g input.gfa -o output -f node
````

### align 


### mod
Modify a graph or alignment based plink file (bed, bim, fam). Use gretl to export nodes or paths which should be excluded from the ped file. 
