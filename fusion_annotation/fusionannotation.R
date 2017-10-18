############################################################
# see examples of Rscript for executing fusionannotation_exec.R
# Rscript --vanilla fusionannotation_exec.R starfusion_output /home/exacloud/lustre1/users/patterja/Functions/git fusions_annotated.txt
# 
############################################################
require(data.table)

starfusion_to_bedpe <- function(starfusion){
  #' starfusion output to bedpe
  #' @description
  #' convert star-fusion.fusion_candidates.final.abridged.FFPM to a bedpe
  #' @param starfusion star-fusion.fusion_candidates.final.abridged.FFPM
  #' @return data.table in bedpe format
  #' @examples
  #' starfusion="/.../star-fusion.fusion_candidates.final.abridged.FFPM"
  #' fusion_bedpe_dt=starfusion_to_bedpe(starfusion)
  #' @export
  starfusion = fread(starfusion)
  left=do.call('rbind',strsplit(as.character(starfusion$LeftBreakpoint), split=":"))
  right=do.call('rbind',strsplit(as.character(starfusion$RightBreakpoint), split=":"))
  leftgene=do.call('rbind',strsplit(as.character(starfusion$LeftGene), split="\\^"))
  rightgene=do.call('rbind',strsplit(as.character(starfusion$RightGene), split="\\^"))
  
  fusion_bedpe_dt=cbind(left[,1],as.numeric(left[,2])-1, as.numeric(left[,2]), right[,1], as.numeric(right[,2])-1, as.numeric(right[,2]), paste(leftgene[,1],"-", rightgene[,1], sep=""),0,left[,3],right[,3],starfusion[,c(2:4,14,15,9)], leftgene[,2], starfusion[,c(10:11)], rightgene[,2], starfusion[,c(12:13)],leftgene[,1],rightgene[,1])
  colnames(fusion_bedpe_dt)=c("chrom1","start1","end1","chrom2","start2","end2","name","score","strand1","strand2","JunctionReadCount","SpanningFragCount","SpliceType","J_FFPM","S_FFPM","LargeAnchorSupport","LeftGeneEns","LeftBreakDinuc","LeftBreakEntropy","RightGeneEns","RightBreakDinuc","RightBreakEntropy","leftgene","rightgene")
  
  return(fusion_bedpe_dt)
  
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
annoOncotator <- function(){
  #' @param
  #' @param
  
  annoOncotator
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
combineOncoKBrefs <- function(annovarsPath=NULL, actionvarsPath=NULL){
  #' combine ONCOKB references
  #' @description
  #' Combine an OncoKB annotation files (allAnnotatedVariants.txt and allActionableVariants.txt) to a combined annotation file, returns a combined data.table of ONCOKB variants.
  #' @param annovarsPath string path to allAnnotatedVariants.txt
  #' @param actionvarsPath string path to allActionableVariants.txt
  #' @return data table of combined ONCOKB variants \code{annovarsPath} and \code{actionvarsPath}
  #' @examples
  #' actionvarsPath=paste0(resource_vardir,"OncoKB/030117/allActionableVariants.txt")
  #' annovarsPath=paste0(resource_vardir,"OncoKB/030117/allAnnotatedVariants.txt")
  #' allvars_oncokb_dt= combineOncoKBrefs(annovarsPath = annovarsPath, actionvarsPath = actionvarsPath)
  #' @export
  # import oncokb tables
  annovars_dt=fread(file=annovarsPath, quote="")
  actionvars_dt=fread(file=actionvarsPath)
  actionvars_dt$table="actionable"
  colnames(actionvars_dt)=c("Gene", "Alteration", "Cancer.Type.or.Oncogenicity", "Level.or.MutationEffect", "Drugs.s.","PMIDs.for.drug.or.mutationeffect","Abstracts.for.drug.or.mutationeffect","table")
  #make annovar columns like actionvars
  annovars_dt$Drugs.s="NA"
  annovars_dt$table="annotated"
  #rename
  colnames(annovars_dt)=c("Gene", "Alteration", "Cancer.Type.or.Oncogenicity", "Level.or.MutationEffect", "PMIDs.for.drug.or.mutationeffect","Abstracts.for.drug.or.mutationeffect","Drugs.s.","table")
  #reorder
  setcolorder(annovars_dt,c("Gene", "Alteration", "Cancer.Type.or.Oncogenicity", "Level.or.MutationEffect", "Drugs.s.","PMIDs.for.drug.or.mutationeffect","Abstracts.for.drug.or.mutationeffect","table"))
  
  allvars_dt=rbind(annovars_dt, actionvars_dt)
  return(allvars_dt)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Annotate with OncoKB
annoOncoKB_fusions <- function(allvars_oncokb_dt, fusion_bedpe_dt){
  #' Annotate a bedpe with OncoKB combined allvars_oncokb_dt
  #' 
  #' @param allvars_oncokb_dt data.table combined ONCOKB variants with combineOncoKBrefs function
  #' @param starfusion fusion_bedpe_dt data.table bedpe format with starfusion_to_bedpe function
  #' @return fusion_bedpe_annoOncokb, data.table annotated with oncokb, cbind to right of bedpe
  #' overlap is based on matching Hugo_Symbol? Not so good but all we have on Oncokb
  #' WANT TO FIX THIS ANNOTATION ONLY inverted and exact matches, no left/right only matches
  #' @examples
  #' fusion_bedpe_dt=starfusion_to_bedpe(starfusion)
  #' actionvarsPath=paste0(resource_vardir,"OncoKB/030117/allActionableVariants.txt")
  #' annovarsPath=paste0(resource_vardir,"OncoKB/030117/allAnnotatedVariants.txt")
  #' allvars_oncokb_dt= combineOncoKBrefs(annovarsPath = annovarsPath, actionvarsPath = actionvarsPath)
  #' fusion_bedpe_annoOncokb_dt=annoOncoKB_fusions(allvars_oncokb_dt = allvars_oncokb_dt, fusion_bedpe_dt = fusion_bedpe_dt)
  #' @export
  # FORMAT ONKOKB:
  #1) get only fusion events
  fusion_oncokb_dt=allvars_oncokb_dt[grep(".*Fusion*",allvars_oncokb_dt$Alteration),]
  #2) get combined fusion names where they exist and replace the single gene name if combo name exists
  #example: combo name= "GPIAP1-PDGFRB"
  alteration_comboindex=grep(" Fusion", fusion_oncokb_dt$Alteration)
  alteration_comboname=do.call('rbind',strsplit(as.character(fusion_oncokb_dt$Alteration[alteration_comboindex]), split=" "))[,1]
  fusion_oncokb_dt$Gene[alteration_comboindex]=alteration_comboname
  #3) aggregate to unique fusion events so no duplicate annotation is necessary
  #aggregation to string
  fusion_oncokb_uniq_dt=fusion_oncokb_dt[ , lapply(.SD, toString), by = list(Gene)]
  
  # create empty data.table
  #the data=NA prevents warnings with data type later (string)
  oncokbannotated_dt=data.table(matrix(nrow=nrow(fusion_bedpe_dt), ncol=ncol(fusion_oncokb_uniq_dt), data="NA")) 
  colnames(oncokbannotated_dt)=paste("ONCOKB", colnames(fusion_oncokb_uniq_dt), sep = "_")
  
  
  for (i in 1:nrow(fusion_bedpe_dt)) {
    if (fusion_bedpe_dt$name[i] %in% fusion_oncokb_uniq_dt$Gene){
      print(1)
      print(paste("index",i))
      annorow=as.list(fusion_oncokb_uniq_dt[fusion_oncokb_uniq_dt$Gene==fusion_bedpe_dt$name[i],])
      oncokbannotated_dt[i, names(oncokbannotated_dt) := annorow]
    }
    else if (fusion_bedpe_dt$leftgene[i] %in% fusion_oncokb_uniq_dt$Gene) {
      print(2)
      print(paste("index",i))
      annorow=as.list(fusion_oncokb_uniq_dt[fusion_oncokb_uniq_dt$Gene==fusion_bedpe_dt$leftgene[i],])
      oncokbannotated_dt[i, names(oncokbannotated_dt) := annorow]
    }
    else if (fusion_bedpe_dt$rightgene[i] %in% fusion_oncokb_uniq_dt$Gene) {
      print(3)
      print(paste("index",i))
      annorow=as.list(fusion_oncokb_uniq_dt[fusion_oncokb_uniq_dt$Gene==fusion_bedpe_dt$rightgene[i],])
      oncokbannotated_dt[i, names(oncokbannotated_dt) := annorow]
    }
  }
  fusion_bedpe_annoOncokb=cbind(fusion_bedpe_dt, oncokbannotated_dt)
  return(fusion_bedpe_annoOncokb)
}


annoCivic_fusions <- function(fusion_bedpe_dt, resource_vardir){
  #' Annotate Civic database
  #' @description Civic table is auto uploaded 01-Mar-2017-ClinicalEvidenceSummaries.tsv
  #' @param starfusion fusion_bedpe_dt data.table bedpe format 
  #' @param resource_vardir In BioCoders resource_vardir="/home/users/patterja/BioCoders/DataResources/Variation/"
  #' Here it is fusion_bedpe_dt (from starfusion_to_bedpe) or fusion_bedpe_annoOncokb_dt (from combineOncoKBrefs)
  #'
  #' @return fusion_bedpe_annoCivic, data.table of annotated with civic, cbind to right of bedpe
  #' overlap is based on overlap of genomic coordinates and matching Hugo_Symbol
  #' @examples
  #' fusion_bedpe_dt=starfusion_to_bedpe(starfusion)
  #' fusion_bedpe_annoCivic_dt=annoCivic_fusions(fusion_bedpe_dt=fusion_bedpe_annoOncokb_dt, resource_vardir)
  #' @export
  # IMPORT CIVIC
  civic=fread(input = paste(resource_vardir,"Civic/030117/01-Mar-2017-ClinicalEvidenceSummaries.tsv", sep=""))
  #civic=fread("/Users/patterja/Workspace/FusionAnnotation/Civic/030117/01-Mar-2017-ClinicalEvidenceSummaries.tsv")
  
  #Some erros in Civic where end is smaller than the start.
  mistake=civic[438,.(stop)]
  civic[438,stop:=(civic[438,.(start)])]
  civic[438,start:=mistake]
  mistake=civic[128,.(stop2)]
  civic[128,stop2:=(civic[128,.(start2)])]
  civic[128,start2:=mistake]

  #get fusions only
  fusion_civic_dt=civic[grep("AMPLIFICATION", variant, invert=TRUE)][!is.na(start2)]

  
  #get unique combo-gene names only
  combo_name_ls=sapply(strsplit(fusion_civic_dt$variant, split=" "),`[`,1)
  dups_ls=names(table(combo_name_ls))[table(combo_name_ls)>1] # replicates
  
  #replace variant names with replicates with "Type1" and weird annotation
  for (i in dups_ls){
    fusion_civic_dt$variant[grep(paste(i,"*", sep=""),fusion_civic_dt$variant)]=i
  }
  #AGGREGATION by variant (there are duplicates for contradicting evidence and weird duplicates)
  lsUniq_func <-{function(x) toString(unique(x))} 
  fusion_civic_byvar_dt=fusion_civic_dt[ , lapply(.SD, lsUniq_func), by = list(variant)]
  
  #some start and stop overlapping duplicates
  fusion_civic_byvar_dt$start=do.call(rbind, lapply(fusion_civic_byvar_dt$start, min))[,1]
  fusion_civic_byvar_dt$start2=do.call(rbind, lapply(fusion_civic_byvar_dt$start2, min))[,1]
  fusion_civic_byvar_dt$stop=do.call(rbind, lapply(fusion_civic_byvar_dt$stop, max))[,1]
  fusion_civic_byvar_dt$stop2=do.call(rbind, lapply(fusion_civic_byvar_dt$stop2, max))[,1]
  
  # Used to do by GenomicRanges overlaps, see older version in gitlab. 
  # Now overlap of same gene name ONLY
  
  # create empty data.frame
  civicannotated_df= data.frame(matrix(nrow=nrow(fusion_bedpe_dt), ncol=ncol(fusion_civic_byvar_dt))) #the data=NA prevents warnings later
  colnames(civicannotated_df)=paste("CIVIC",colnames(fusion_civic_byvar_dt), sep="_")

  
  for (i in 1:nrow(fusion_bedpe_dt)) {
    if (fusion_bedpe_dt$name[i] %in% fusion_civic_byvar_dt$variant){
      civicannotated_df[i,]=data.frame(fusion_civic_byvar_dt[fusion_civic_byvar_dt$variant==fusion_bedpe_dt$name[i],])
    }
    else if (paste0(fusion_bedpe_dt$rightgene[i],"-",fusion_bedpe_dt$leftgene[i]) %in% fusion_civic_byvar_dt$variant){
      civicannotated_df[i,]=data.frame(fusion_civic_byvar_dt[fusion_civic_byvar_dt$variant==paste0(fusion_bedpe_dt$rightgene[i],"-",fusion_bedpe_dt$leftgene[i]),])
    } else {
      civicannotated_df[i,]=rep("NA", ncol(fusion_civic_byvar_dt))
    }
  }
  
  
  fusion_bedpe_annoCivic=cbind(data.table(fusion_bedpe_dt), data.table(civicannotated_df))
  return(fusion_bedpe_annoCivic)
}

  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  