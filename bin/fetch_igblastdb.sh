#!/usr/bin/env bash
# Download IgBLAST database files
#
# Author:  Jason Vander Heiden
# Date:    2020.08.11
#
# Arguments:
#   -o = Output directory for downloaded files. Defaults to current directory.
#   -x = Include download of internal_data and optional_file bundles.
#   -h = Display help.

# Default argument values
OUTDIR="."
DOWNLOAD_ALL=false

# Print usage
usage () {
    echo "Usage: `basename $0` [OPTIONS]"
    echo "  -o  Output directory for downloaded files. Defaults to current directory."
    echo "  -x  Include download of legacy internal_data and optional_file bundles."
    echo "  -h  This message."
}

# Get commandline arguments
while getopts "o:xh" OPT; do
    case "$OPT" in
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
        ;;
    x)  DOWNLOAD_ALL=true
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

# Make output directory if it does not exist
if $OUTDIR_SET && [ ! -d "${OUTDIR}" ]; then
    mkdir -p $OUTDIR
fi

# Fetch database
wget -q -r -nH --cut-dirs=5 --no-parent \
    ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/database \
    -P ${OUTDIR}/database

# Extract
tar -C ${OUTDIR}/database -xf ${OUTDIR}/database/mouse_gl_VDJ.tar
tar -C ${OUTDIR}/database -xf ${OUTDIR}/database/rhesus_monkey_VJ.tar

if $DOWNLOAD_ALL; then
    # Fetch internal_data
    wget -q -r -nH --cut-dirs=5 --no-parent \
        ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/old_internal_data \
        -P ${OUTDIR}/internal_data

    # Fetch optional_file
    wget -q -r -nH --cut-dirs=5 --no-parent \
        ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/old_optional_file  \
        -P ${OUTDIR}/optional_file
fi