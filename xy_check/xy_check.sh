#!/bin/bash

VERSION='2.0.0'

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
    echo "{\"xy_check\": 0}" > "output.txt";
    XY="FEMALE"
elif [ $COUNT -gt $3 ]; then
    echo "{\"xy_check\": 1}" > "output.txt";
    XY="MALE"
else
    echo "{\"xy_check\": 999}" > "output.txt";
    XY="UNSPECIFIED"
fi
echo $XY >> "log.txt"
