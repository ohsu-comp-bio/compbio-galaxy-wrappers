#!/bin/bash

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
    echo "Gender does not match Samplesheet"
    echo "Gender does not match Samplesheet" >&2
    exit 1
fi
