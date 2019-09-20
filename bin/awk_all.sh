#!/bin/bash
PREFIX="clonal_analysis/Clone_lineage/Graphml_trees/noheader" # store outputs
FOLDER="clonal_analysis/Clone_lineage/Graphml_trees"   # input files
mkdir -p "${PREFIX}"    # make sure the outputs dir exists

for FILE in "${FOLDER}/*.txt"       # get the file names you want to work on
do
  # use ${PREFIX}/${FILE} to redirect output to a 
  # file that's associated with the input
  awk '/<graph /,/<\/graph>/' "${FILE}" > "${PREFIX}/${basename -- $FILE}"
done

