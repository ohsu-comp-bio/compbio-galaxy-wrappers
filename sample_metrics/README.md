
# Sample Metrics v0.9.2   

**USAGE**  
sample_metrics.py [-h]  
                 [--probeqc_after PROBEQC_AFTER]  
                 [--probeqc_before PROBEQC_BEFORE]  
                 [--fastqc_r1 FASTQC_R1]  
                 [--fastqc_r2 FASTQC_R2]  
                 [--picard_summary PICARD_SUMMARY]  
                 [--picard_summary_umi PICARD_SUMMARY_UMI]  
                 [--gatk_depth_cov_prop GATK_DEPTH_COV_PROP]  
                 [--gatk_depth_cov_cnts GATK_DEPTH_COV_CNTS]  
                 [--gatk_coll_rnaseq_mets GATK_COLL_RNASEQ_METS]  
                 [--gatk_count_reads_total GATK_COUNT_READS_TOTAL]  
                 [--gatk_count_reads_ints GATK_COUNT_READS_INTS]  
                 [--msi MSI]  
                 [--primers_bam PRIMERS_BAM]  
                 [--primers_bed PRIMERS_BED]  
                 [--blia_pre BLIA_PRE]  
                 [--blia_post BLIA_POST]  
                 [--dragen_metrics DRAGEN_METRICS]  
                 [--dragen_qc DRAGEN_QC]  
                 [--json_in [JSON_IN [JSON_IN ...]]]  
                 [--outfile OUTFILE]  
                 [--outfile_txt OUTFILE_TXT]  
                 [--version]  

optional arguments:  
  -h, --help            show this help message and exit  
  --probeqc_after PROBEQC_AFTER
                        Probe coverage QC after UMI deduplication metrics.  
  --probeqc_before PROBEQC_BEFORE
                        Probe coverage QC before UMI deduplication metrics.  
  --fastqc_r1 FASTQC_R1
                        FastQC stats for read 1.  
  --fastqc_r2 FASTQC_R2
                        FastQC stats for read 2.  
  --picard_summary PICARD_SUMMARY
                        Picard alignment summary metrics file.  
  --picard_summary_umi PICARD_SUMMARY_UMI
                        Picard alignment summary metrics file taken after umi
                        deduplication.  
  --gatk_depth_cov_prop GATK_DEPTH_COV_PROP
                        GATK DepthOfCoverage file.  
  --gatk_depth_cov_cnts GATK_DEPTH_COV_CNTS
                        GATK DepthOfCoverage locus file.  
  --gatk_coll_rnaseq_mets GATK_COLL_RNASEQ_METS
                        GATK CollectRnaSeqMetrics file.  
  --gatk_count_reads_total GATK_COUNT_READS_TOTAL
                        Output from GATK4 CountReads with total read count.  
  --gatk_count_reads_ints GATK_COUNT_READS_INTS
                        Output from GATK4 CountReads with read count for reads
                        overlapping targets.  
  --msi MSI             TSV file containing MSI results  
  --primers_bam PRIMERS_BAM
                        BAM file to calculate primer reads on target.  
  --primers_bed PRIMERS_BED
                        BED file containing primer coordinates only.  
  --blia_pre BLIA_PRE   JSON from Ding correlation subtyping, pre-
                        normalization.  
  --blia_post BLIA_POST
                        JSON from Ding correlation subtyping, post-
                        normalization.  
  --dragen_metrics DRAGEN_METRICS
                        JSON file containing metrics produced by DRAGEN  
  --dragen_qc DRAGEN_QC
                        CSV file containing fastqc data produced by DRAGEN  
  --json_in [JSON_IN [JSON_IN ...]]
                        Arbitrary number of files to be included in sample
                        metrics that are in json format.  
  --outfile OUTFILE     Output file with json string.  
  --outfile_txt OUTFILE_TXT
                        Output file in human readable text format.
  --version             show program's version number and exit  


----
**METRICS**  
qthirty: The percentage of bases with a quality score of 30 or higher. Q30 means 1 in 1000 or 99.9% probability of incorrect base call  
qthirty_before: qthirty before deduplication 
averageDepth: 
averageDepth_before  
depthTen  
depthTen_before  
depthTwenty  
depthTwenty_before  
depthFifty  
depthFifty_before  
depthOneHundred  
depthOneHundred_before  
depthTwoHundredFifty  
depthTwoHundredFifty_before  
depthFiveHundred  
depthFiveHundred_before  
depthSevenHundred  
depthSevenHundred_before  
depthOneThousand  
depthOneThousand_before  
depthTwelveHundredFifty  
depthTwelveHundredFifty_before  
depthTwoThousand  
depthTwoThousand_before  
uniformity_of_coverage: Percentage of sites with coverage greater than 20% of the mean coverage in region.  
pf_reads  
pf_reads_after  
pct_pf_reads_aligned:  
pct_pf_reads_aligned_after  
allele_balance  
allele_balance_het_count  
gatk_cr_on_target  
gatk_cr_total  
gatk_cr_ints  
gatk_pct_mrna_bases  
gatk_pct_correct_strand_reads  
gc_pct_r1  
gc_pct_r2  
bio_sex_check  
homozygosity_flag  
parentage_binom  
parentage_disc  
parentage_confirmed  
parentage_sites  
percentOnTarget  
percentOnTarget_after  
percentUmi  
rna_count_zero  
rna_count_one  
rna_count_ten  
rna_count_onehundred  
rna_count_onethousand  
rna_count_tenthousand  
rna_count_hundredthousand  
rna_tpm_zero  
rna_tpm_hundredth  
rna_tpm_tenth  
rna_tpm_one  
rna_tpm_ten  
rna_tpm_onehundred  
rna_tpm_onethousand  
tmb  
msi_pct  
msi_sites  
msi_somatic_sites  
blia_pre_best  
blia_pre_second  
blia_pre_diff  
blia_post_best  
blia_post_second  
blia_post_diff  
blia_reportable  
blia_raw_pre  
blis_raw_pre  
lar_raw_pre  
mes_raw_pre  
blia_raw_post  
blis_raw_post  
lar_raw_post  
mes_raw_post  
total_on_target_transcripts  
total_on_target_transcripts_pct  
forced_calls_above  
forced_calls_below  
y_ploidy_check  
cnv_median_segment_mad_cn  
q30_bases_pct  
average_alignment_coverage_over_target_region  
pct_of_target_region_with_coverage_10x_inf  
pct_of_target_region_with_coverage_20x_inf  
pct_of_target_region_with_coverage_50x_inf  
pct_of_target_region_with_coverage_100x_inf  
aligned_reads_in_target_region_pct  
dragen_gc_pct_r1  
dragen_gc_pct_r2  
number_of_large_roh_gt_eq_3000000  
ploidy_estimation  


----
**VERSION HISTORY**   
0.9.2  
- Update to pass list of metrics to be provided. (Collected from config file during workflow generation.)  

0.9.1  
- Update to allow for case when picard_summary and picard_summary_umi are None  

0.9.0  
- Remove old style metrics in favor of new style  

0.8.15  
- Change IlluminaExome_V2_PlusMito to IlluminaExome_V2-5_PlusMito  

0.8.14  
- Update tmb to use 0.0 as placeholder value  
- Update QIAseq_V4_STP4 to remove bio_sex_check metric  

0.8.13  
- Add pre-umi dedup depth metrics and post-umi alignment metrics  
- Changed workflow name to IlluminaExome_V2_PlusMito  

0.8.12  
- Add class DragenMetrics and class DragenQC to handle metrics data for IlluminaExome_V2_DRAGEN  

0.8.11  
- Add QIAseq_XP_RNA_STP  

0.8.10  
- Switched before and after picard metric source in _pumi  

0.8.9  
- Add uniformity_of_coverage, cnv_median_segment_mad_cn metrics  

0.8.8  
- Add QIAseq_V4_STP4  
- Revert xy_check to bio_sex_check to maintain CGD compatibility  

0.8.7  
- Rename metric bio_sex_check to y_ploidy_check for AgilentCRE_V1 and TruSightOneV2_5  
- Rename metric bio_sex_check to xy_check for QIAseq_V4_MINI  

0.8.5  
- Change forced calls metric to json  

0.8.3  
- Add entry for json metric bio_sex_check  

0.8.0  
- Pass forced calls above and below background metric  

0.7.0  
- Now calculate percentOnTarget from FastQC tatal sequences and collectalignmentmetrics on sorted bwa bam  


