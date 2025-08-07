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
REPERTOIRE="imgt"
DATE=$(date +"%Y.%m.%d")
TMPDIR=$(mktemp -d)

# Associative array (for BASH v3) where keys are species folder names and values are query strings
# TODO: get back to human
#SPECIES_QUERY=("human:Homo%20sapiens"
#                "mouse:Mus%20musculus")
SPECIES_QUERY=("mouse:Mus%20musculus")
# Associative array (for BASH v3) with species name replacements
SPECIES_REPLACE=('human:s/Homo sapiens/Homo_sapiens/g'
                'mouse:s/Mus musculus/Mus_musculus/g')

# Counter for loop iteration, used for getting the right values of SPECIES_REPLACE
COUNT=0
# For each species
for SPECIES in ${SPECIES_QUERY[@]}
do
    if [ "$SPECIES" == "human:Homo%20sapiens" ]; then
        STRAINS=("")
    elif [ "$SPECIES" == "mouse:Mus%20musculus" ]; then
        STRAINS=("C57BL-6:C57BL%25252f6")
    fi

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
        if [ "$KEY" == "human" ]; then
            echo "|------ ${LOCUS}"
            URL="https://ogrdb.airr-community.org/download_germline_set/${VALUE}/${LOCUS}_${GENES}/published/gapped_ex"
            FILE_NAME="${TMPDIR}/${REPERTOIRE}_${KEY}_${LOCUS}.fasta"
            wget $URL -O $FILE_NAME -q
        elif [ "$KEY" == "mouse" ]; then
            for STRAIN in ${STRAINS[@]}
            do
                STRAIN_KEY=${STRAIN%%:*}
                STRAIN_VALUE=${STRAIN#*:}
                echo "|------- ${STRAIN_KEY}"
                echo "|-------- ${LOCUS}"
                # Adding ugly case as C57BL-6 AIRR-C reference does not exist for IGK and IGL and using C57BL-6J instead
                if [ "$LOCUS" == "IGH" ]; then
                    URL="https://ogrdb.airr-community.org/download_germline_set/${VALUE}/${STRAIN_VALUE}/${STRAIN_VALUE}%20${LOCUS}/published/gapped"
                    echo $URL
                    FILE_NAME="${TMPDIR}/${REPERTOIRE}_${KEY}_${LOCUS}.fasta"
                    wget $URL -O $FILE_NAME -q
                elif [ "$LOCUS" == "IGK" ] || [ "$LOCUS" == "IGL" ]; then
                    # C57BL-6J is the strain used in AIRR-C for IGK and IGL
                    # C57BL-6 is not available in AIRR-C for IGK and IGL
                    if [ "$STRAIN_KEY" == "C57BL-6" ]; then
                        STRAIN_VALUE="C57BL%25252f6J"
                    fi

                    # Get V gene data and check it downloaded well
                    URLV="https://ogrdb.airr-community.org/download_germline_set/${VALUE}/${STRAIN_VALUE}/${STRAIN_VALUE}%20${LOCUS}V/published/gapped"
                    echo $URLV
                    FILENAME_V="${FILE_PATH}/${REPERTOIRE}_${KEY}_${LOCUS}V.fasta"
                    wget $URLV -O $FILENAME_V -q
                    if [ -s "$FILENAME_V" ]
                    then
                        # Check that file starts with ">"
                        if ! grep -q "^>" "$FILENAME_V"; then
                            echo "Error: Fasta file ${FILENAME_V} does not start with '>'"
                            cat ${FILENAME_V} | head -n 10
                            exit 1
                        fi
                    else
                        echo "AIRR-C Fasta file does not exist, or is empty. Is the OGRDB server online?"
                        exit 1
                    fi

                    # Get J gene data and check it downloaded well
                    URLJ="https://ogrdb.airr-community.org/download_germline_set/${VALUE}/${LOCUS}J%20(all%20strains)/published/gapped"
                    echo $URLJ
                    FILENAME_J="${FILE_PATH}/${REPERTOIRE}_${KEY}_${LOCUS}J.fasta"
                    wget $URLJ -O $FILENAME_J -q
                    if [ -s "$FILENAME_J" ]
                    then
                        # Check that file starts with ">"
                        if ! grep -q "^>" "$FILENAME_J"; then
                            echo "Error: Fasta file ${FILENAME_J} does not start with '>'"
                            cat ${FILENAME_J} | head -n 10
                            exit 1
                        fi
                    else
                        echo "AIRR-C Fasta file does not exist, or is empty. Is the OGRDB server online?"
                        exit 1
                    fi

                else
                    STRAIN_VALUE=${STRAIN#*:}
                fi

                #https://ogrdb.airr-community.org/download_germline_set/Mus%20musculus/C57BL%25252f6/C57BL%25252f6%20IGH/published/gapped
                #TODO: change how to handle mouse strains in folder structure
                #TODO: right now it will only download the files for the first strain and treat them as "mouse"
                #KEY="${KEY}_${STRAIN_KEY}"
            done
        fi



        for SEGMENT in ${SEGMENTS[@]}
        do
            if [ "$KEY" == "mouse" ] && ([ "$LOCUS" == "IGK" ] || [ "$LOCUS" == "IGL" ]); then
                continue
            fi
            # Split the file into segments
            F="${REPERTOIRE}_${KEY}_${LOCUS}${SEGMENT}.fasta"
            FILE_NAME="${FILE_PATH}/${F}"
            ./split_airrcdb.py ${TMPDIR}/${REPERTOIRE}_${KEY}_${LOCUS}.fasta ${FILE_NAME} ${LOCUS} ${SEGMENT}

            # Checking once that file exists and is not empty (checks IMGT server is online)
            if [ -s "$FILE_NAME" ]
            then
                # Check that file starts with ">"
                if ! grep -q "^>" "$FILE_NAME"; then
                    echo "Error: Fasta file ${FILE_NAME} does not start with '>'"
                    cat ${FILE_NAME} | head -n 10
                    exit 1
                fi
            else
                echo "AIRR-C Fasta file does not exist, or is empty. Is the OGRDB server online?"
                exit 1
            fi
        done


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

    # Set back to the original KEY for mouse
    if [ "$SPECIES" == "mouse:Mus%20musculus" ]; then
        KEY=${SPECIES%%:*}
    fi
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
