#!/bin/bash

"""
Check the heterozygosity of X to confirm specification in samplesheet

author: JL
"""

gatk IndexFeatureFile -F $1

gatk SelectVariants -V $1 -O selected.vcf --select-type-to-include SNP --exclude-filtered true -L X -sample-name $2 -select "vc.getGenotype(\"$2\").isHet()"

gatk VariantsToTable -V selected.vcf -O output_HOP.tsv -F CHROM -F POS -GF GT

read lines filename <<< $(wc -l output_HOP.tsv)

if [ $lines -lt 2 ]; then
    BIO_SEX="Male"
    echo "Male"
else
    BIO_SEX="Female"
    echo "Female"
fi

if [ $BIO_SEX == $3 ]; then
    echo "Biological sex matches." >> output.txt
else
    COUNT=$[lines-1]
    echo "Biological sex does NOT match. Samplesheet specifies $3, but there are $COUNT heterozygous spots on X, indicating that the sample is likely $BIO_SEX." >> output.txt
    echo "Biological sex does NOT match. Samplesheet specifies $3, but there are $COUNT heterozygous spots on X, indicating that the sample is likely $BIO_SEX." >&2
    
    exit 1
fi