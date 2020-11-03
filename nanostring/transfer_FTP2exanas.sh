#!/bin/bash
## Usage:
## transfers ftp to Xdrive ../HSR/CLINICAL/Nanostring/automated_data directory

wf_version="1.0.0"
REPO_DIR="/home/users/${USER}/nanostring"

config=${REPO_DIR}"/nanostring_config.json"
STORAGE_HOST=$( jq -r  '.exanas.host' "${config}" )
STORAGE_GROUP=$( jq -r  '.exanas.group' "${config}" )
TODAY=$(date +"%Y%m%d")

echo $config
echo $STORAGE_HOST
echo $STORAGE_GROUP

# X DIRECTORY SYNC
#cd "${X_DATA_DIR}"

# tabs and spaces seem to mess with heredocs
## Get the data using heredoc
#ftp -inv $RAW_HOST > backupRCC_"${TODAY}".log 2>&1 <<EOF
#user $RAW_USER $RAW_PASSWORD
#cd RCCData
#binary
#mget --no-clobber *
#echo "mgetting"
#ls
#bye
#EOF
echo "running"
ssh "${STORAGE_HOST}" "sg ${STORAGE_GROUP} \"echo "running script" && python ${REPO_DIR}"/"transfer_FTP2exanas.py --config ${REPO_DIR}"/"nanostring_config.json 2>&1 | tee /home/users/patterja/log.file\""


