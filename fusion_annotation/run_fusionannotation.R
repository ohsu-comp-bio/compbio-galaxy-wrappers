#!/usr/bin/env Rscript
############################################################
# HOW TO RUN
# Rscript --vanilla run_fusionannotation.R star-fusion.fusion_candidates.final.abridged.FFPM oncotated_left_output oncotated_right_output left_sequence.txt right_sequence.txt
# 6/26/2016 Adjusted for galaxy input flags
############################################################
#location of DataResources
#resource_vardir="/home/exacloud/lustre1/BioCoders/DataResources/Variation/"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ParseArguments: probably not the best way to do this with positional arguments. Use argparse() argparser() but this works for now as a galaxy wrapper
args = commandArgs(trailingOnly=TRUE)

starfusion = args[1]
oncotated_left_output = args[2]
oncotated_right_output = args[3]
sequence_left = args[4]
sequence_right = args[5]
annovarsPath = args[6]
actionvarsPath = args[7]
civicTablePath = args[8]

#starfusion="/Volumes/exacloud/BioCoders/ProjectCollaborations/SMMART_RNA_Workflows/fusion_detection/output_mpssr/RNA161026CC_KDL-RNA-Fusion-1_S1/star-fusion.fusion_candidates.final.abridged.FFPM"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (length(args)==0) {
  stop("Must supply arguments! 
       Usage: Rscript --vanilla run_fusionannotation.R star-fusion.fusion_candidates.final.abridged.FFPM \
       oncotated_left_output \
       oncotated_right_output \
       left_sequence.txt \
       right_sequence.txt \
       /home/exacloud/lustre1/BioCoders/DataResources/Variation/OncoKB/030117/allAnnotatedVariants.txt \
       /home/exacloud/lustre1/BioCoders/DataResources/Variation/OncoKB/030117/allActionableVariants.txt", call.=FALSE)
} else if (length(args)==5) {
    annovarsPath=paste0("allAnnotatedVariants.txt")
    actionvarsPath=paste0("allActionableVariants.txt")
    civicTablePath=paste0("01-Mar-2017-ClinicalEvidenceSummaries.tsv")
  }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Load Libraries

library(data.table)
source("/home/groups/clinical/installedTest/galaxy-dist/tools/fusion_annotation/fusionannotation.R")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Prep anno files
#see combineOncoKBrefs function

allvars_oncokb_dt= combineOncoKBrefs(annovarsPath = annovarsPath, actionvarsPath = actionvarsPath)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fusion_bedpe_dt=starfusion_to_bedpe(starfusion)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#add oncotator output

#direct_name=paste0(dirname(starfusion), "/anno_files/")
oncotated_left=fread(skip = 1,oncotated_left_output)
oncotated_right=fread(skip = 1,oncotated_right_output) 

#select columns
onco_addleft=oncotated_left[,c("HGNC_RefSeq IDs","transcript_exon", "COSMIC_FusionGenes_fusion_genes"), with=FALSE]
onco_addright=oncotated_right[,c("HGNC_RefSeq IDs","transcript_exon", "COSMIC_FusionGenes_fusion_genes"), with=FALSE]
names(onco_addleft)=paste("oncotatorleft", names(onco_addleft), sep = "_")
names(onco_addright)=paste("oncotatorright", names(onco_addright), sep = "_")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#add sequence output
seqleft=fread(sequence_left ,header=FALSE)
seqright=fread(sequence_right ,header=FALSE) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fusion_bedpe_dt_oncotatseq=cbind(fusion_bedpe_dt, onco_addleft, onco_addright, seq_left=seqleft$V2, seq_right=seqright$V2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#annotate with OncoKB, see annoOncoKB_fusions function
fusion_bedpe_annoOncokb_dt=annoOncoKB_fusions(allvars_oncokb_dt = allvars_oncokb_dt, fusion_bedpe_dt = fusion_bedpe_dt_oncotatseq)

fusion_bedpe_annoCivic_dt=annoCivic_fusions(fusion_bedpe_dt=fusion_bedpe_annoOncokb_dt, civic_table = civicTablePath)

#write.table(, file=args[2], row.names=FALSE)
write.table(fusion_bedpe_annoCivic_dt, file="fusions_annotated.txt", row.names=FALSE, quote=FALSE, sep="\t")
  
