#!/bin/bash

#curl -X GET -d "user_token=iSUQvGmVNSjq834g9fP5" -d "user_email=campbena@ohsu.edu" https://oregon.moon.diploid.com/samples/$1/patient-info > info.txt
python moon_api.py -id $1 -mode info

HPOS=$( grep 'showterm?id=HP:' info.txt | sed "s/[<][^>]*[>]//g")

COUNT=1

for LINE in $HPOS; do
  if [[ $LINE == *HP:* ]]; then
    if [[ $COUNT == 1 ]]; then
      COUNT=2
      echo -n "$LINE" > hpo.txt
    else
      echo -n ";$LINE" >> hpo.txt
    fi
  fi
done

SELECTEDS=$(grep 'selected="selected"' info.txt)

GENDER="unknown"

for LINE in $SELECTEDS; do
  if [[ $LINE == *male* ]]; then
    GENDER="male"
  fi
  if [[ $LINE == *female* ]]; then
    GENDER="female"
  fi
done

echo $GENDER

AGE=$(grep 'id="sample_age"' info.txt | cut -d '"' -f 8)

echo $AGE

CONSANG=$(grep 'id="sample_is_consanguinous"' info.txt | cut -d '"' -f 6)

IS_CONSANG="false"

if [[ $CONSANG == 1 ]]; then
  let IS_CONSANG="true"
fi

echo $IS_CONSANG

TRIO_CHECK=$(grep 'value="mother"' info.txt)

if [[ ! -z $TRIO_CHECK ]]; then
  TRIO=$(grep 'selected="selected"' info.txt)
  ID="nothing"
  REL="nothing"
  HEALTH="nothing"
  NUM=1
  LAST_LINE="nothing"
  TWO_LAST="nothing"
  IFS=$'\n'
  for LINE in $TRIO; do
    echo "The line is:" $LINE
      THREE_LAST=$TWO_LAST
      echo $THREE_LAST
      TWO_LAST=$LAST_LINE
      echo $TWO_LAST
      LAST_LINE=$LINE
      echo $LAST
    if [[ $LINE == *KD* ]]; then
      ID=$(echo $LINE | cut -d '"' -f 4 )
      echo $ID
    elif [[ $LINE == *healthy* ]]; then
      HEALTH="healthy"
      echo $HEALTH
    elif [[ $LINE == *affected* ]]; then
      HEALTH="affected"
      echo $HEALTH
    elif [[ $LINE == *mother* ]]; then
      REL="mother"
      echo $REL
      ID=$(echo $THREE_LAST | cut -d '"' -f 4 )
      if [[ $NUM == 2 ]]; then
        echo -n ";$ID:$HEALTH:$REL" >> parents.txt
      else
        echo -n "$ID:$HEALTH:$REL" > parents.txt
      fi
      let NUM=2
    elif [[ $LINE == *father* ]]; then
      REL="father"
      echo $REL
      ID=$(echo $THREE_LAST | cut -d '"' -f 4 )
      if [[ $NUM == 2 ]]; then
        echo -n ";$ID:$HEALTH:$REL" >> parents.txt
      else
        echo -n "$ID:$HEALTH:$REL" > parents.txt
      fi
      let NUM=2
    fi
  done
fi

HPO_TERMS=$(less hpo.txt)

PARENTS=$(less parents.txt)

NEW_ID=$1

if [[ ! -z $TRIO_CHECK ]]; then
echo "  curl -F ?user_token=iSUQvGmVNSjq834g9fP5?\ "
echo "       -F ?user_email=campbena@ohsu.edu?\ "
echo "       -F ?snp_vcf_file=@$2?\ "
echo "       -F ?sv_vcf_file=@$3?\ "
echo "       -F ?age=$AGE?\ "
echo "       -F ?gender=$GENDER?\ "
echo "       -F ?is_consanguinous=$IS_CONSANG?\ "
echo "       -F ?hpo_terms=$HPO_TERMS?\ "
echo "       -F ?family_members=$PARENTS?\ "
echo "       https://oregon.moon.diploid.com/samples.json"

python moon_api.py -snp $2 -cnv $3 -mode post -age $AGE -gender $GENDER -consang $IS_CONSANG -hp $HPO_TERMS -family $PARENTS
  #NEW_ID=$(curl -F "user_token=j4yKzsZivJuCXzoxdUM6"\
       #-F "user_email=potteram@ohsu.edu"\
       #-F "snp_vcf_file=@$2"\
       #-F "sv_vcf_file=@$3"\
       #-F "age=$AGE"\
       #-F "gender=$GENDER"\
       #-F "is_consanguinous=$IS_CONSANG"\
       #-F "hpo_terms=$HPO_TERMS"\
       #-F "family_members=$PARENTS"\
       #https://oregon.moon.diploid.com/samples.json | sed 's/[^0-9]*//g' )
else

  python moon_api.py -snp $2 -cnv $3 -mode post -age $AGE -gender $GENDER -consang $IS_CONSANG -hp $HPO_TERMS
  
#  NEW_ID=$(curl -F "user_token=j4yKzsZivJuCXzoxdUM6"\
#       -F "user_email=potteram@ohsu.edu"\
#       -F "snp_vcf_file=@$2"\
#       -F "sv_vcf_file=@$3"\
#       -F "age=$AGE"\
#       -F "gender=$GENDER"\
#       -F "is_consanguinous=$IS_CONSANG"\
#       -F "hpo_terms=$HPO_TERMS"\
#       https://oregon.moon.diploid.com/samples.json | sed 's/[^0-9]*//g' )

fi
 
NEW_ID=$(sed 's/[^0-9]*//g' new_id.txt)
echo $NEW_ID

python moon_api.py -id $NEW_ID -mode analyse

#curl -F "user_token=iSUQvGmVNSjq834g9fP5" -F "user_email=campbena@ohsu.edu" https://oregon.moon.diploid.com/samples/$NEW_ID/analysis.json       
       
       
       