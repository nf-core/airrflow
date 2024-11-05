#!/usr/bin/env bash
# set -euo pipefail
# Convert AIRR-C and IMGT germline sequences to IgBLAST database
#
# Author:  Gisela Gabernet
# Date:    2024.11.03
# Licence: AGPL-3
#
# Arguments:
#   -i = Input directory containing germlines in the form <species>/vdj/imgt_<species>_<LOCUS><segment>.fasta
#   -o = Output directory for the built database. Defaults to current directory.
#   -h = Display help.

# Default argument values

OUTDIR="."

# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -i  Input directory containing germlines in the form:"
    echo -e "      <species>/vdj/airrc-imgt_<species>_<chain><segment>.fasta."
    echo -e "  -o  Output directory for the built database."
    echo -e "  -h  This message."
}



# Get commandline arguments
while getopts "i:o:h" OPT; do
    case "$OPT" in
    i)  GERMDIR=$(realpath $OPTARG)
        GERMDIR_SET=true
        ;;
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
        ;;
    h)  usage
        exit
        ;;
    \?) echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    :)  echo "Option -$OPTARG requires an argument" >&2
        exit 1
        ;;
    esac
done

# Exit if no germline directory provided
if ! $GERMDIR_SET; then
    echo "You must specify an input directory using the -i option" >&2
    exit 1
fi

# Create and set directories
OUTDIR=$(realpath ${OUTDIR})
mkdir -p ${OUTDIR}/fasta
TMPDIR=$(mktemp -d)

# Create fasta files of each species, LOCUS and segment combination
for SPECIES in human #mouse
do
    for CHAIN in IG TR
    do
        for SEGMENT in V D J
        do
            F=$(echo airrc-imgt_${SPECIES}_${CHAIN}_${SEGMENT}.fasta | tr '[:upper:]' '[:lower:]')
            cat ${GERMDIR}/${SPECIES}/vdj/airrc-imgt_${SPECIES}_${CHAIN}?${SEGMENT}.fasta > ${TMPDIR}/${F}
        done

        # C nucleotides
        F=$(echo airrc-imgt_${SPECIES}_${CHAIN}_c.fasta | tr '[:upper:]' '[:lower:]')
        cat ${GERMDIR}/${SPECIES}/constant/airrc-imgt_${SPECIES}_${CHAIN}?C.fasta > ${TMPDIR}/${F}
    done
done

# Parse each created fasta file to create igblast database
cd ${TMPDIR}
NT_FILES=$(ls *.fasta | grep -E "airrc-imgt_(human|mouse)_ig_(v|d|j)\.fasta")
echo ${NT_FILES}
for F in ${NT_FILES}; do
    cp ${F} ${OUTDIR}/fasta/${F}
    makeblastdb -parse_seqids -dbtype nucl -in ${OUTDIR}/fasta/${F} \
        -out ${OUTDIR}/database/${F%%.*}
done

# Reference data from IMGT needs cleaning of the headers
C_FILES=$(ls *.fasta | grep -E "airrc-imgt_(human|mouse)_(ig|tr)_c\.fasta")
TR_FILES=$(ls *.fasta | grep -E "airrc-imgt_(human|mouse)_tr_(v|d|j)\.fasta")
for F in ${IMGT_FILES[@]}; do
    clean_imgtdb.py ${F} ${OUTDIR}/fasta/${F}
    makeblastdb -parse_seqids -dbtype nucl -in ${OUTDIR}/fasta/${F} \
        -out ${OUTDIR}/database/${F%%.*}
done

# Remove temporary fasta files
cd -; rm -rf $TMPDIR
