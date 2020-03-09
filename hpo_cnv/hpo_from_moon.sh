#!/bin/bash

echo "email:$2"

echo "token:$3"

HPOS=$(curl -X GET -d "user_token=$3" -d "user_email=$2" https://oregon.moon.diploid.com/samples/$1/patient-info | grep 'showterm?id=HP:' | sed "s/[<][^>]*[>]//g")

COUNT=1

for LINE in $HPOS; do
  if [[ $LINE == *HP:* ]]; then
    if [[ $COUNT == 1 ]]; then
      COUNT=2
      echo "$LINE" > hpo.txt
    else
      echo "$LINE" >> hpo.txt
    fi
  fi
done