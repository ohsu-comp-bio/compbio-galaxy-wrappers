#!/bin/bash

mkdir $2
cr=$'\r'

while read -r ID; do

ID="${ID%$cr}"

PRE='https://hpo.jax.org/api/hpo/download/term?identifier='

POST='&association=genes'

URL="${PRE}${ID}${POST}"

curl -XGET $URL --output $2/$ID.xlsx

done < $1