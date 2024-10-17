#!/bin/bash
###############################################################################
# Verify that the graph has paths, not walks, and that the path names
# confirm to the Pan SN spec, and that the 
#       [sample_name][delim][haplotype_id][delim][contig_or_scaffold_name]
# where [delim] is a delimiter character, typically a #.
#
# Usage: check_graph.bash <graph.gfa> <delim>
###############################################################################

# Enable bash strict mode
set -uo pipefail

# Check that the correct number of arguments were provided
if [[ "$#" -ne 2 ]]; then
    echo "Usage: check_graph.bash <graph.gfa> <delim>"
    exit 1
fi

# Extract the arguments
GRAPH_GFA=$1
DELIM=$2


# Check that there are paths at all
if [[ "$(grep '^P' $GRAPH_GFA | head -1 )" == "" ]] ; then
    echo "Error: the graph ${GRAPH_GFA} does not have any paths."
    exit 1
fi
echo "✓ Graph has paths."

# Check that the graph has paths, not walks
if [[ "$(grep '^W' $GRAPH_GFA |head -1)" != "" ]] ; then
    echo "Error: The graph ${GRAPH_GFA} has walks. These should be converted to paths first."
    echo "TODO: Script/info about converting walks to paths (vg something)."
    exit 1
fi
echo "✓ Graph has no walks."


# Check that the path names are unique -- the path names are in the 2nd column
# of the GFA file
if [[ "$(awk '$1 == "P" {print $2}' $GRAPH_GFA | sort | uniq -c | awk '$1 > 1' | grep . |head -1 )" != "" ]] ; then
    echo "Error: The path names are not unique."
    exit 1
fi
echo "✓ Path names are unique."


# Check that the path names confirm to the Pan SN spec
# [sample_name][delim][haplotype_id][delim][contig_or_scaffold_name]
# where
# sample_name := string
# delim := character
# haplotype_id := number
# contig_or_scaffold_name := string
# The path names are in the 2nd column of the GFA file
INCORRECT_PATH_NAMES="$(awk '$1 == "P" {print $2}' $GRAPH_GFA | grep -Ev '^.+'$DELIM'[0-9]+'$DELIM'.+$' | grep . |head -1 )"
if [[ "${INCORRECT_PATH_NAMES}" != "" ]] ; then
    echo "Error: The path names do not confirm to the Pan SN spec."
    echo "First incorrect path name: ${INCORRECT_PATH_NAMES}"
    echo "The path names should confirm to the Pan SN spec:"
    echo "  [sample_name][delim][haplotype_id][delim][contig_or_scaffold_name]"
    echo "Current delimiter: ${DELIM}"
    echo "More information: https://github.com/pangenome/PanSN-spec"
    exit 1
fi
echo "✓ Path names confirm to the Pan SN spec."