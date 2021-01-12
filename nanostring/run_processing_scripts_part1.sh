#/bin/bash
#
# Pipeline to run everything together. 
## Usage as loop
## for batch in /Volumes/OHSU/CLINICAL/Nanostring/automated_data/*; do echo $batch; run_processing_scripts.sh $batch; done
##
## Usage:
## run_processing_scripts.sh 20190314_208421591020
##

valid_version="20200512"
ihc_version="20200507"
knownpos_version="1.0"
antibody_version="1.0"

#second flag is true/false will include or exclude blanks 
NEWBATCH=$(basename "$1")
omit_blank=$2

echo $NEWBATCH
DATA_DIR="/Volumes/Histopathology\ Shared\ Resource/CLINICAL/Nanostring"
REPO_DIR="/Users/letaw/PycharmProjects/compbio-galaxy-wrappers/nanostring"

ss_dir="${DATA_DIR}""/""automated_data/""$NEWBATCH"
output_dir="${DATA_DIR}""/""Assay_Whole_Slide/output/""$NEWBATCH"
#files
FILES=("${DATA_DIR}""/""automated_data/"$NEWBATCH"/*RCC")
abfile="${DATA_DIR}""/""Assay_Whole_Slide/REFERENCE_FILES/ANTIBODY_REFERENCE_v"${antibody_version}".csv"
known_pos_file="${DATA_DIR}""/""Assay_Whole_Slide/REFERENCE_FILES/knownpositives_v"${knownpos_version}".txt"
validationraw_file="${DATA_DIR}""/""Assay_Whole_Slide/REFERENCE_FILES/validation_samples_rawdata_"${valid_version}".txt"
ihc_file="${DATA_DIR}""/""Assay_Whole_Slide/REFERENCE_FILES/ihc_status_"${ihc_version}".txt"
md_file="/Users/patterja/Box/NANOSTRING/nanostring_metadata.xlsx"

#find samplesheet
echo $NEWBATCH
echo "${FILES[@]}"
ss=$(find "${ss_dir}" -name "${NEWBATCH}*samplesheet.txt")
echo $ss

#Step 1: if samplesheet exists process raw RCC files
if test -f "$ss"; then
    mkdir -p "${output_dir}" 
    echo "dir made"
    cd "${output_dir}" 
    echo `pwd`
    if [ "$omit_blank" = "true" ]; then
        echo "omitting blank RCCs"
        /usr/local/bin/python3 "${REPO_DIR}""/"process_rcc.py --omit_blank --samplesheet "$ss" --abfile "${abfile}" ${FILES[@]}
    else
        /usr/local/bin/python3 "${REPO_DIR}""/"process_rcc.py --samplesheet "$ss" --abfile "${abfile}" ${FILES[@]}
    fi
else
    echo "no samplesheet"
    exit 1
fi

#Step 2: if rawdata.txt exists then check that qc_metrics are within range
if [ -s rawdata.txt ]; then
    /usr/local/bin/python3 "${REPO_DIR}""/"check_qc.py --rawdata "rawdata.txt" --runmetrics "run_metrics.txt"
    echo "finished checking qc"
else
    echo "something went wrong with processing rawdata.txt or run_metrics.txt"
    exit 1
fi

#Step 3: Check qc flags in json
echo "check qc flags"
echo $PWD
if [ -e qc_metrics.json ]; then
    echo "entering"
    pf=$(cat qc_metrics.json | jq '.[] | .[] | .qc' | sort | uniq -c | grep "FAIL")
    if [[ ! -z "$pf" ]]; then #check if any pass/fail flag is not empty, something failed
        echo $NEWBATCH
        echo "check qc_metrics.json. QC flag failed"
        if [ "$omit_blank" = "true" ]; then
           echo "keepting rawdata.txt as is. CHECK QC FLAGS AND SAMPLES!"
        else
          mv rawdata.txt rawdata_QCFLAG.txt
        fi
    else
        echo "All QC PASS"
    fi
else
    echo "check that qc_metrics.json is formed correctly"
    exit 1
fi

#Step 4: if rawdata.txt exists process QC stats for batch pass/fail
if [ -s rawdata.txt ]; then
    "${REPO_DIR}""/"process_batch_qc.R -i "rawdata.txt" --pos_file "${known_pos_file}" --validation_file "${validationraw_file}" --md_file "${md_file}" 
else
    echo "something went wrong with making rawdata.txt"
    exit 1
fi

echo "Process 1 complete. Check qc_positive_controls.tsv for PASS/FAIL flags."

