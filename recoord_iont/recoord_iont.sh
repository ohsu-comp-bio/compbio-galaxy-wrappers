#!/bin/bash

# Modify VCF coordinates produced by TorrentSuite to play nice with Richards lab workflow reference genome hg19.
# recoord_iont.sh <old_file.vcf> <new_file.vcf> <old_file.bam> <new_file.sam>

SAMTOOLS="/opt/installed/samtools"

cat $1 | sed 's/chr//g' | sed -e 's/^M/MT/g' | sed 's/ID=M,/ID=MT,/g' > $2

$SAMTOOLS view -h $3 | sed 's/chrM/chrMT/g' | sed 's/chr//g' > $4
#$SAMTOOLS view -b $4 > $5
