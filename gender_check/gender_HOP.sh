#!/bin/bash

gatk IndexFeatureFile -F $1

gatk SelectVariants -V $1 -O selected.vcf --select-type-to-include SNP --exclude-filtered true -L X -sample-name $2 -select "vc.getGenotype(\"$2\").isHet()"

gatk VariantsToTable -V selected.vcf -O output_HOP.tsv -F CHROM -F POS -GF GT

read lines filename <<< $(wc -l output_HOP.tsv)

if [ $lines -lt 2 ]; then
    GENDER="Male"
    echo "Male"
else
    GENDER="Female"
    echo "Female"
fi

if [ $GENDER == $3 ]; then
    echo "Gender matches." >> output.txt
else
    COUNT=$[lines-1]
    echo "Gender does NOT match. Samplesheet has gender as $3, but there are $COUNT heterozygous spots on X, indicating that the sample is likely $GENDER." >> output.txt
    echo "Gender does NOT match. Samplesheet has gender as $3, but there are $COUNT heterozygous spots on X, indicating that the sample is likely $GENDER." >&2
    
    exit 1
fi