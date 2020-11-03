#!/bin/bash
#
# Pipeline to RUV batch correction . 
## Usage as loop
## for batch in /Volumes/OHSU/CLINICAL/Nanostring/automated_data/*; do echo
##$batch; run_processing_scripts_part2.sh $batch; done
## --no-cohort: exclude comparison with MBC
##
## Usage:
## run_processing_scripts_part2.sh 20190314_208421591020
##
set -e

valid_version="20200512"
ihc_version="20200507"
knownpos_version="1.0"
antibody_version="1.0"

#cd "/Users/patterja/Box\ Sync/NANOSTRING/data/${smmart}"
NEWBATCH=$(basename "$1")
echo $NEWBATCH
DATA_DIR="/Volumes/OHSU/CLINICAL/Nanostring"
REPO_DIR="/Users/patterja/Workspace/nanostring/nanostring"

ss_dir="${DATA_DIR}""/""automated_data/""$NEWBATCH"
output_dir="${DATA_DIR}""/""Assay_Whole_Slide/output/""$NEWBATCH"
#files
FILES=("${DATA_DIR}""/""automated_data/"$NEWBATCH"/*RCC")
abfile="${DATA_DIR}""/""Assay_Whole_Slide/REFERENCE_FILES/ANTIBODY_REFERENCE_v"${antibody_version}".csv"
known_pos_file="${DATA_DIR}""/""Assay_Whole_Slide/REFERENCE_FILES/knownpositives_v"${knownpos_version}".txt"
validationraw_file="${DATA_DIR}""/""Assay_Whole_Slide/REFERENCE_FILES/validation_samples_rawdata_"${valid_version}".txt"
ihc_file="${DATA_DIR}""/""Assay_Whole_Slide/REFERENCE_FILES/ihc_status_"${ihc_version}".txt"
md_file="/Users/patterja/Box/NANOSTRING/nanostring_metadata.xlsx"

#Step 4: check output directory exists and step1 was run 

if [ ! -d "${output_dir}" ]; then
    echo "something went wrong. Check batch id and output directory location"
    exit 1
fi
# change directory to output dir
cd "${output_dir}"
echo $PWD

# Step 5:  if rawdata.txt exists process sample batch correction with cohort  and batch QC has been reviewed
# Step 6: write report for each sample
if [ -s rawdata.txt ]; then
    echo "running batch correction"
#    /Users/patterja/Workspace/nanostring/nanostring/batch_correction_ruv.R -i "rawdata.txt" --validation_file "${validationraw_file}" --md_file "${md_file}" --ab_ref_file "${abfile}"
    "${REPO_DIR}""/"batch_correction_ruv.R --input "rawdata.txt" --validation_file "${validationraw_file}" --md_file "${md_file}" --ihc_file "${ihc_file}" --ab_ref_file "${abfile}"
#    run report
    "${REPO_DIR}""/"nanostring_report.R --repo "${REPO_DIR}" --rawdata "rawdata.txt" --md_file "${md_file}"


else
    echo "something went wrong with rawdata file,\n make sure you ran run_processing_script_part1.sh"
    exit 1
fi
