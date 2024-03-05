# gfa2bin

Convert a gfa to a plink file. This tool can also use a compressed coverage file from packing.  

## Subcommands
### graph - Converting from graphs directly

You are able to convert a gfa directly to a plink file using the inherent graph structure. This does include nodes, edges, and directed nodes. We support PanSN-spec and it can help to easily merge path the genome sample. 

#### Example usage: 
````text
gfa2bin graph -g input.gfa -o output -f node
````

### align 


### mod
Modify a graph or alignment based plink file (bed, bim, fam). Use gretl to export nodes or paths which should be excluded from the ped file. 
