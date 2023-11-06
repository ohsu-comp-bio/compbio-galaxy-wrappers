#!/bin/bash

VERSION='2.1.0'

gatk CollectReadCounts -I $1 -L Y -O output.tsv -imr OVERLAPPING_ONLY -format TSV

COUNT=$(
python << END

with open("output.tsv", "r") as count_file:
    count = 0
    for line in count_file:
        if line.startswith("Y"):
            line_array = line.split("\t")
            count = int(line_array[3])
            print(count)
END
)
echo $COUNT > "log.txt"

if [ $COUNT -lt $2 ]; then
    echo "{\"bio_sex_check\": \"Female\"}" > "output.txt";
    GENDER="FEMALE"
elif [ $COUNT -gt $3 ]; then
    echo "{\"bio_sex_check\": \"Male\"}" > "output.txt";
    GENDER="MALE"
else
    echo "{\"bio_sex_check\": \"Indeterminate\"}" > "output.txt";
    GENDER="UNSPECIFIED"
fi
echo $GENDER >> "log.txt"
