#!/bin/bash

DIR="Graphml_trees"   # input files

IDX=0
for FILE in ${DIR}/*.graphml       # get the file names you want to work on
do
    if [ $IDX -eq 0 ]
    then
        awk '/<graph id=/ {exit} {print}' $FILE > "${DIR}/head.txt"
        tail -n 1 $FILE > "${DIR}/tail.txt"
        awk '/<graph /,/<\/graph>/' "${FILE}" >> "${DIR}/All_graphs_patient_nohead.txt"
        ((IDX++))
    else
        awk '/<graph /,/<\/graph>/' "${FILE}" >> "${DIR}/All_graphs_patient_nohead.txt"
        ((IDX++))
    fi
done
cat "${DIR}/head.txt" "${DIR}/All_graphs_patient_nohead.txt" "${DIR}/tail.txt" > "${DIR}/All_graphs_patient.graphml"
rm "${DIR}/All_graphs_patient_nohead.txt" "${DIR}/head.txt" "${DIR}/tail.txt"

