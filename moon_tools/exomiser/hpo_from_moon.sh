#!/bin/bash


HPOS=$(curl -X GET -d "user_token=$2" -d "user_email=$1" https://oregon.moon.diploid.com/samples/$1/patient-info | grep 'showterm?id=HP:' | sed "s/[<][^>]*[>]//g")

COUNT=1

for LINE in $HPOS; do
  if [[ $LINE == *HP:* ]]; then
    if [[ $COUNT == 1 ]]; then
      COUNT=2
      echo -n "'$LINE'" > hpo.txt
    else
      echo -n ",'$LINE'" >> hpo.txt
    fi
  fi
done
