#!/usr/bin/env bash
# Written by Gisela Gabernet and released under the MIT license (2020).
echo "Fetching databases..."

fetch_imgt.sh -o imgtdb_base

fetch_igblastdb.sh -x -o igblast_base

imgt2igblast.sh -i ./imgtdb_base -o igblast_base

echo "FetchDBs process finished."
