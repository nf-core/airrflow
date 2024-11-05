#!/usr/bin/env bash
# Download germlines from the AIRR-C and IMGT website
# So far the AIRR-C reference provides only IG VDJ regions
# This script downloads IG VDJ regions from AIRR-C and IG constant regions from IMGT
# As well as TCR VDJ and constant regions from IMGT
#
# Author:  Mohamed Uduman, Jason Anthony Vander Heiden, Gisela Gabernet
# Date:    2024.11.04
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
REPERTOIRE="airrc-imgt"
DATE=$(date +"%Y.%m.%d")
TMPDIR=$(mktemp -d)

# Associative array (for BASH v3) where keys are species folder names and values are query strings
SPECIES_QUERY=("human:Homo%20sapiens")
               # "mouse:Mus")
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
    mkdir -p $FILE_PATH

    # VDJ Ig
    echo "|---- Ig"
    for LOCUS in IGH IGK IGL
    do
        if [ "$LOCUS" == "IGH" ]; then
            GENES="VDJ"
            SEGMENTS=("V" "D" "J")
        else
            GENES="VJ"
            SEGMENTS=("V" "J")
        fi
        URL="https://ogrdb.airr-community.org/download_germline_set/${VALUE}/${LOCUS}_${GENES}/published/gapped_ex"
        FILE_NAME="${TMPDIR}/${REPERTOIRE}_${KEY}_${LOCUS}.fasta"
        wget $URL -O $FILE_NAME -q

        for SEGMENT in ${SEGMENTS[@]}
        do
            F="${REPERTOIRE}_${KEY}_${LOCUS}${SEGMENT}.fasta"
            split_airrcdb.py ${TMPDIR}/${REPERTOIRE}_${KEY}_${LOCUS}.fasta ${FILE_PATH}/${F} ${LOCUS} ${SEGMENT}
        done

        # Checking once that file exists and is not empty (checks IMGT server is online)
        if [ -s "$FILE_NAME" ]
        then
            echo "AIRR-C Fasta file exists and is not empty"
        else
            echo "AIRR-C Fasta file does not exist, or is empty. Is the OGRDB server online?"
            exit 1
        fi
    done

    # VDJ TCR
    echo "|---- TCR"
    for CHAIN in TRAV TRAJ TRBV TRBD TRBJ TRDV TRDD TRDJ TRGV TRGJ
    do
        URL="https://www.imgt.org/genedb/GENElect?query=7.14+${CHAIN}&species=${VALUE}"
        FILE_NAME="${FILE_PATH}/${REPERTOIRE}_${KEY}_${CHAIN}.fasta"
        TMP_FILE="${FILE_NAME}.tmp"
        #echo $URL
        wget $URL -O $TMP_FILE -q
        awk '/<pre>/{i++}/<\/pre>/{j++}{if(j==2){exit}}{if(i==2 && j==1 && $0!~"^<pre>"){print}}' $TMP_FILE > $FILE_NAME
        sed -i.bak "$REPLACE_VALUE" $FILE_NAME && rm $FILE_NAME.bak
        rm $TMP_FILE
    done

    # Download constant regions from IMGT (AIRR-C does not provide C region references)
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

        URL="https://www.imgt.org/genedb/GENElect?query=${QUERY}+${CHAIN}&species=${VALUE}"
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
        URL="https://www.imgt.org/genedb/GENElect?query=14.1+${CHAIN}&species=${VALUE}"
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

# Write download info for AIRR-C
INFO_FILE=${OUTDIR}/IMGT.yaml
echo -e "source:  https://www.imgt.org/genedb" > $INFO_FILE
echo -e "date:    ${DATE}" >> $INFO_FILE
echo -e "species:" >> $INFO_FILE
for Q in ${SPECIES_QUERY[@]}
do
    echo -e "    - ${Q}" >> $INFO_FILE
done

# Write download info for IMGT (Cregions)
INFO_FILE=${OUTDIR}/AIRRC.yaml
echo -e "source:  https://ogrdb.airr-community.org/download_germline_set" > $INFO_FILE
echo -e "date:    ${DATE}" >> $INFO_FILE
echo -e "species:" >> $INFO_FILE
for Q in ${SPECIES_QUERY[@]}
do
    echo -e "    - ${Q}" >> $INFO_FILE
done
