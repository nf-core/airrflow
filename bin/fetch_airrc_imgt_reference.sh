#!/usr/bin/env bash
# set -euo pipefail
# Written by Gisela Gabernet and released under the MIT license (2020).

echo "Fetching databases..."

./fetch_airrc_imgt.sh -o airrc_imgt_base

./fetch_igblastdb.sh -x -o igblast_base

./airrc_imgt2igblast.sh -i ./airrc_imgt_base -o igblast_base

echo "FetchDBs process finished."
