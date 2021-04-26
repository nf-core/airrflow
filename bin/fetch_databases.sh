    echo "Fetching databases..."

    fetch_imgt.sh -o imgtdb_base

    fetch_igblastdb.sh -x -o igblast_base

    imgt2igblast.sh -i ./imgtdb_base -o igblast_base

    echo "FetchDBs process finished."