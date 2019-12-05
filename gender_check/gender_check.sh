#!/bin/bash

VERSION='1.1.1'

gatk CollectReadCounts -I $2 -L Y -O output.tsv -imr OVERLAPPING_ONLY -format TSV

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
echo $COUNT
echo $COUNT > "output.txt"

if [ $COUNT -lt $3 ]; then
    GENDER="FEMALE"
elif [ $COUNT -gt $4 ]; then
    GENDER="MALE"
else
    GENDER="UNSPECIFIED"
fi
echo $GENDER
echo $GENDER >> "output.txt"

if [ $GENDER = $1 ]; then
    echo "Gender matches Sample sheet"
else
    if [ -z ${5} ]; then
        echo "Gender does not match Samplesheet"
    else
        cat output.txt | mailx -s "$(hostname) Gender Check Tool Error" "${5}"
    fi
    echo "Gender does not match Samplesheet" >&2
    exit 1
fi
