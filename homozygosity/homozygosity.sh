#!/bin/bash

gatk IndexFeatureFile -F $1

gatk SelectVariants -select "QUAL>250" -V $1 -O selected.vcf --select-type-to-include SNP --exclude-filtered true

gatk VariantsToTable -V selected.vcf -O output.tsv -F CHROM -F POS -GF GT

Rscript gtrellis_homozygosity.R output.tsv

python3 zygosity_runs.py output.tsv > nums.txt