#!/usr/bin/env bash
# Written by Gisela Gabernet and released under the MIT license (2020).
set -euo pipefail

DATABASE_TYPE=imgt

while getopts "d:" OPT; do
    case "$OPT" in
    d)  DATABASE_TYPE=$OPTARG
        ;;
    esac
done

echo "Fetching databases..."

bash fetch_references.sh -d "${DATABASE_TYPE}" -o reference_base

fetch_igblastdb.sh -x -o igblast_base

bash ref2igblast.sh -d "${DATABASE_TYPE}" -i ./reference_base -o igblast_base

echo "FetchDBs process finished."
