#!/usr/bin/env bash
# Download germlines from the IMGT website
#
# Author:  Mohamed Uduman, Jason Anthony Vander Heiden
# Date:    2017.07.03
# Licence: AGPL-3
#
# Arguments:
#   -o = Output directory for downloaded files. Defaults to current directory.
#   -h = Display help.

# Default argument values
OUTDIR="."

# Print usage
usage () {
    echo "Usage: `basename $0` [OPTIONS]"
    echo "  -o  Output directory for downloaded files. Defaults to current directory."
    echo "  -h  This message."
}

# Get commandline arguments
while getopts "o:h" OPT; do
    case "$OPT" in
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

# Info
REPERTOIRE="imgt"
DATE=$(date +"%Y.%m.%d")

# Associative array (for BASH v3) where keys are species folder names and values are query strings
SPECIES_QUERY=("human:Homo+sapiens"
                "mouse:Mus")
# Associative array (for BASH v3) with species name replacements
SPECIES_REPLACE=('human:s/Homo sapiens/Homo_sapiens/g'
                'mouse:s/Mus musculus/Mus_musculus/g')

# Counter for loop iteration, used for getting the right values of SPECIES_REPLACE
COUNT=0
# For each species
for SPECIES in ${SPECIES_QUERY[@]}
do
    KEY=${SPECIES%%:*}
    VALUE=${SPECIES#*:}
    REPLACE_VALUE=${SPECIES_REPLACE[$COUNT]#*:}
    echo "Downloading ${KEY} repertoires into ${OUTDIR}"

    # Download VDJ
    echo "|- VDJ regions"
    FILE_PATH="${OUTDIR}/${KEY}/vdj"
    FILE_PATH_AA="${OUTDIR}/${KEY}/vdj_aa"
    mkdir -p $FILE_PATH $FILE_PATH_AA

    # VDJ Ig
    echo "|---- Ig"
    for CHAIN in IGHV IGHD IGHJ IGKV IGKJ IGLV IGLJ
    do
        URL="http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.14+${CHAIN}&species=${VALUE}"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME

        # Checking once that file exists and is not empty (checks IMGT server is online)
        read file
        if [ -s "$FILE_NAME" ]
        then
            echo "IMGT Fasta file exists and is not empty"
        else
            echo "IMGT Fasta file does not exist, or is empty. Is the IMGT server online?"
            exit 1
        fi

        # Make sed command work also for mac, see: https://stackoverflow.com/a/44864004
        sed -i.bak "$REPLACE_VALUE" $FILE_NAME && rm $FILE_NAME.bak
        rm $TMP_FILE
    done

    # V amino acid for Ig
    for CHAIN in IGHV IGKV IGLV
    do
        URL="http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+${CHAIN}&species=${VALUE}"
        FILE_NAME="${FILE_PATH_AA}/${REPERTOIRE}_aa_${KEY}_${CHAIN}.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        # Make sed command work also for mac, see: https://stackoverflow.com/a/44864004
        sed -i.bak "$REPLACE_VALUE" $FILE_NAME && rm $FILE_NAME.bak
        rm $TMP_FILE
    done

    # VDJ TCR
    echo "|---- TCR"
    for CHAIN in TRAV TRAJ TRBV TRBD TRBJ TRDV TRDD TRDJ TRGV TRGJ
    do
        URL="http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.14+${CHAIN}&species=${VALUE}"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        sed -i.bak "$REPLACE_VALUE" $FILE_NAME && rm $FILE_NAME.bak
        rm $TMP_FILE
    done

    # V amino acid for TCR
    for CHAIN in TRAV TRBV TRDV TRGV
    do
        URL="http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+${CHAIN}&species=${VALUE}"
        FILE_NAME="${FILE_PATH_AA}/${REPERTOIRE}_aa_${KEY}_${CHAIN}.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        sed -i.bak "$REPLACE_VALUE" $FILE_NAME && rm $FILE_NAME.bak
        rm $TMP_FILE
    done

    # Download leaders
    echo "|- Spliced leader regions"
    FILE_PATH="${OUTDIR}/${KEY}/leader"
    mkdir -p $FILE_PATH

    # Leader Ig
    echo "|---- Ig"
    for CHAIN in IGH IGK IGL
    do
        URL="http://www.imgt.org/IMGT_GENE-DB/GENElect?query=8.1+${CHAIN}V&species=${VALUE}&IMGTlabel=L-PART1+L-PART2"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}L.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        sed -i.bak "$REPLACE_VALUE" $FILE_NAME && rm $FILE_NAME.bak
        rm $TMP_FILE
    done

    # Leader TCR
    echo "|---- TCR"
    for CHAIN in TRA TRB TRG TRD
    do
        URL="http://www.imgt.org/IMGT_GENE-DB/GENElect?query=8.1+${CHAIN}V&species=${VALUE}&IMGTlabel=L-PART1+L-PART2"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}L.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        sed -i.bak "$REPLACE_VALUE" $FILE_NAME && rm $FILE_NAME.bak
        rm $TMP_FILE
    done

    # Download constant regions
    echo "|- Spliced constant regions"
    FILE_PATH="${OUTDIR}/${KEY}/constant/"
    mkdir -p $FILE_PATH

    # Constant Ig
    echo "|---- Ig"
    for CHAIN in IGHC IGKC IGLC
    do
        QUERY=14.1
        if [ "${KEY}" == "mouse" ] && ([ "$CHAIN" == "IGKC" ] || [ "$CHAIN" == "IGLC" ]); then
            # IMGT does not have artificially spliced mouse IGKC / IGLC
            QUERY=7.5
        fi

        URL="http://www.imgt.org/IMGT_GENE-DB/GENElect?query=${QUERY}+${CHAIN}&species=${VALUE}"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        sed -i.bak "$REPLACE_VALUE" $FILE_NAME && rm $FILE_NAME.bak
        rm $TMP_FILE
    done

    # Constant for TCR
    echo "|---- TCR"
    for CHAIN in TRAC TRBC TRGC TRDC
    do
        URL="http://www.imgt.org/IMGT_GENE-DB/GENElect?query=14.1+${CHAIN}&species=${VALUE}"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        sed -i.bak "$REPLACE_VALUE" $FILE_NAME && rm $FILE_NAME.bak
        rm $TMP_FILE
    done

    echo ""
    ((COUNT++))
done

# Write download info
INFO_FILE=${OUTDIR}/IMGT.yaml
echo -e "source:  http://www.imgt.org/IMGT_GENE-DB" > $INFO_FILE
echo -e "date:    ${DATE}" >> $INFO_FILE
echo -e "species:" >> $INFO_FILE
for Q in ${SPECIES_QUERY[@]}
do
    echo -e "    - ${Q}" >> $INFO_FILE
done
