#!/usr/bin/env bash

# $0 = basename of program!
# $1 = name of run1 file
# $2 = name of run2 file
# $2 = name of run2 file
# $3 = name of  combined trees output file 
# $4 = burninTrees value
# $5 = output MCC tree file

#gets the taxon block from one of the files (in this case the first one) and makes a new file.
grep -v -e "tree STATE" "$1" > "$4"
# for species trees, need to remove the last 'End;'
sed -i '$ d' "$4"
#searches for the trees using the term 'tree STATE' in both tree files and dumps them into alternate lines in the combined file 
paste -d"\n" <(grep "tree STATE" "$1") <(grep "tree STATE" "$2") <(grep "tree STATE" "$3") >> "$4"
#adds the closing nexus code
echo "End;" >> "$4"
# run tree annotator
treeannotator -burninTrees "$5" -heights ca "$4" "$6"

# to run
#./tree_combiner.sh Serra_completo_280717_haps_aligned_RUN1.trees Serra_completo_280717_haps_aligned_RUN2.trees Serra_completo_280717_haps_aligned_RUN3.trees Serra_completo_280717_haps_aligned_COMB.trees 336 Serra_completo_280717_haps_aligned_COMB.tre

# for the tree sample
#treeannotator -burninTrees 0 -heights ca ../Serra_completo_280717_haps_aligned_SAMPLE.trees ../Serra_completo_280717_haps_aligned_SAMPLE.tre
