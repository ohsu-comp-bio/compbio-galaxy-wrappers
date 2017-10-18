#!/bin/bash

if [ "$#" -lt 3 ];
then
  echo "Needs 3 arguments <user/group> <user/group_name> <path>"
  exit -1;
fi

modifier=""
if [ "$1" == "group" ];
then
  modifier="g"
fi

name="$2"
path="$3"

if [ -e "$path" ];
then
  if [ -d "$path" ];
  then
    path="${path}/";
  fi
  #Recursively provide read access
  nfs4_setfacl -RP -a A:${modifier}:${name}@acc.ohsu.edu:rtnc "${path}"
  if [ -d "$path" ];
  then
    #Provide execute access to top level directory
    nfs4_setfacl -P -a A:${modifier}:${name}@acc.ohsu.edu:x "${path}"
    #Provide execute access to all sub-directories
    find ${path} -type d -exec nfs4_setfacl -P -a A:${modifier}:${name}@acc.ohsu.edu:x {} \;
  fi
else
  echo "File/directory \"${path}\" not found"
  exit -1;
fi

echo "Access change complete." > $4
