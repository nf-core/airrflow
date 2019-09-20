#!/usr/bin/bash
PREFIX="noheader"   # define where to store all the outputs
mkdir -p "clonal_analysis/Clone_lineage/Graphml_trees/${PREFIX}"    # make sure the outputs dir exists

for FILE in clonal_analysis/Clone_lineage/Graphml_trees/*.txt       # get the file names you want to work on
do
  # use ${PREFIX}/${FILE} to redirect output to a 
  # file that's associated with the input
  awk '/<graph /,/<\/graph>/' "${FILE}" > "${PREFIX}/${FILE}"
done

