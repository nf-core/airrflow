#!/usr/bin/env bash
# Convert IMGT germlines sequences to IgBLAST database
#
# Author:  Jason Anthony Vander Heiden
# Date:    2016.11.21
#
# Arguments:
#   -i = Input directory containing germlines in the form <species>/vdj/imgt_<species>_<chain><segment>.fasta
#   -o = Output directory for the built database. Defaults to current directory.
#   -h = Display help.

# Default argument values
OUTDIR="."

# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -i  Input directory containing germlines in the form:"
    echo -e "      <species>/vdj/imgt_<species>_<chain><segment>.fasta."
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

# Create fasta files of each species, chain and segment combination
for SPECIES in human mouse
do
    for CHAIN in IG TR
    do
        # VDJ nucleotides
        for SEGMENT in V D J
        do
            F=$(echo imgt_${SPECIES}_${CHAIN}_${SEGMENT}.fasta | tr '[:upper:]' '[:lower:]')
            cat ${GERMDIR}/${SPECIES}/vdj/imgt_${SPECIES}_${CHAIN}?${SEGMENT}.fasta > ${TMPDIR}/${F}
        done

        # V amino acids
        F=$(echo imgt_aa_${SPECIES}_${CHAIN}_v.fasta | tr '[:upper:]' '[:lower:]')
        cat ${GERMDIR}/${SPECIES}/vdj_aa/imgt_aa_${SPECIES}_${CHAIN}?V.fasta > ${TMPDIR}/${F}
    done
done

# Parse each created fasta file to create igblast database
cd ${TMPDIR}
NT_FILES=$(ls *.fasta | grep -E "imgt_(human|mouse).+\.fasta")
for F in ${NT_FILES}; do
    clean_imgtdb.py ${F} ${OUTDIR}/fasta/${F}
    makeblastdb -parse_seqids -dbtype nucl -in ${OUTDIR}/fasta/${F} \
        -out ${OUTDIR}/database/${F%%.*}
done

AA_FILES=$(ls *.fasta | grep -E "imgt_aa_(human|mouse).+\.fasta")
for F in ${AA_FILES}; do
    clean_imgtdb.py ${F} ${OUTDIR}/fasta/${F}
    makeblastdb -parse_seqids -dbtype prot -in ${OUTDIR}/fasta/${F} \
        -out ${OUTDIR}/database/${F%%.*}
done

# Remove temporary fasta files
cd -; rm -rf $TMPDIR
