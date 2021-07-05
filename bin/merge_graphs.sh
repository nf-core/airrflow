#!/bin/bash

DIR="Graphml_trees"   # input files

IDX=0
for FILE in ${DIR}/*.txt       # get the file names you want to work on
do
  if [ $IDX -eq 0 ]
  then
    head -n 26 $FILE > "${DIR}/head.txt"
    tail -n 1 $FILE > "${DIR}/tail.txt"
    awk '/<graph /,/<\/graph>/' "${FILE}" >> "${DIR}/All_graphs_patient_nohead.graphml"
    ((IDX++))
  else
    awk '/<graph /,/<\/graph>/' "${FILE}" >> "${DIR}/All_graphs_patient_nohead.graphml"
    ((IDX++))
  fi
done
cat "${DIR}/head.txt" "${DIR}/All_graphs_patient_nohead.graphml" "${DIR}/tail.txt" > "${DIR}/All_graphs_patient.graphml"
rm "${DIR}/All_graphs_patient_nohead.graphml" "${DIR}/head.txt" "${DIR}/tail.txt"

