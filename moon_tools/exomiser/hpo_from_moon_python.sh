#!/bin/bash


python moon_api.py -mode info -id $1 

HPOS=$(grep 'showterm?id=HP:' info.txt | sed "s/[<][^>]*[>]//g")

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
