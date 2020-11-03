#!/bin/bash
## Usage: transfer_X2exanas.sh automated_data 
## transfers Xdrive to exanas ../HSR/CLINICAL/Nanostring/automated_data directory

wf_version="1.0.0"


config="${1%/}"
#echo $DIR

STORAGE_HOST=$( jq -r  '.exanas.host' "${config}" )
STORAGE_GROUP=$( jq -r  '.exanas.group' "${config}" )
STORAGE_DEST_DIR=$( jq -r  '.exanas.exanas_destdir' "${config}" )"/X"
TODAY=$(date +"%Y%m%d")

X_DATA_DIR=$( jq -r  '.Xdrive.main_dir' "${config}" )
X_SS_DIR==$( jq -r  '.Xdrive.samplesheet_dir' "${config}" )

REMOTE_TEMP="/home/users/patterja/tmp"


echo ${X_DATA_DIR}/ to
echo ${USER}@${STORAGE_HOST}:${REMOTE_TEMP}
 
#move to home directory 
rsync -iv -a --exclude={'archive','logs'} ${X_DATA_DIR}/ ${USER}@${STORAGE_HOST}:${REMOTE_TEMP}
echo files are in remote temp directory

#move from home directory to storage directory
ssh "${STORAGE_HOST}" "sg ${STORAGE_GROUP} \"rsync --remove-source-files -riv ${REMOTE_TEMP}/ ${STORAGE_DEST_DIR}\""

#--remove-source-files doesn't quite work
#print contents of remote tmp directory
ssh ${USER}@${STORAGE_HOST} "ls ${REMOTE_TEMP}"
