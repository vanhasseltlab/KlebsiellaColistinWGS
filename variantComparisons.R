# -	Comparison of all samples against the T=0h sample; 
# o	make a table which lists the number of differential variants found for each of those samples, 
#   and also include the corresponding experimental metadata (experiment type like DCTK, SCTK; time point of sample; replicate ID, strain ID)
# o	an R list object for each sample ID, containing another list with: 
#   - the specific variant identifiers; 
#   - the genes IDs to which the variants are related (if known
#   
# https://pfam.xfam.org/search#tabview=tab1
# https://www.ebi.ac.uk/interpro/

##  LOAD PACKAGES  #########################################################################################################################
pkgs <- c("vcfR", "adegenet", "adegraphics", "pegas", "StAMPP", "lattice", "gplots", "ape", "ggmap", "ggplot2", "reshape","readxl")
pkgs <- c("vcfR", "readxl","writexl","tidyr","stringr","readr")
for (pkg in pkgs){
  if (!pkg %in% installed.packages()) {
    install.packages(pkg, dependencies = TRUE)
  }
}
lapply(pkgs, library, character.only = TRUE)
remove (pkgs,pkg)


##  READ DATA  #########################################################################################################################
# meta data of samples
meta <- read_excel("Sample overview_seq.xlsx")
meta <- as.data.frame(meta)
meta$varCountRaw <- 0
meta$varCountFilt <- 0

strainL <- list()
strainDF <- data.frame()

hyponames <- c("KPN700603" = "NPBKJIEB",
               "KPN9884" = "OLPKOKPN",
               "KPN9749" = "AGIDONFK")

## loop through samples
for (strain in c("KPN9884", "KPN9749", "KPN700603")){
  print (strain)
  strain <- "KPN9749"
  
# read vcf file with vcfR
  vcf <- read.vcfR(paste0("genome_data/toolOut/GATK_pipe_output/",strain,"/annotated.vcf.gz"))

## cleaning ####  
  # remove STK_52 & FSTK_120, sample introduces problems.
  if (strain == "KPN9884"){vcf <- vcf[,head(colnames(vcf@gt),-1)]}
  if (strain == "KPN9749"){vcf <- vcf[,colnames(subset(vcf@gt,select=-c(FSTK_120)))]}
  
  # add varcount per sample to meta
  gt <- extract.gt(vcf, element = "GT",as.numeric = T)
  for (col in colnames(gt)){
    meta[which(meta$`new lable`==col),]$varCountRaw <- sum(gt[,col],na.rm = T)
  }
  
  # get all NG c0 samples and put their genotypes in a dataframe.
  NGsamples <- meta[which(meta$replicate=="NG" & meta$strain==strain & meta$CONC=="c0"),]$`new lable`
  NGvariations <- data.frame(vcf@gt[,c(colnames(vcf@gt)%in%NGsamples)])
  # loop though the columns of the NG dataframe, check if the NG samples have the variation (TRUE/FALSE)
  for (column in NGsamples){
    NGvariations[,column] <- startsWith(vcf@gt[,column], "1:")
  }
  # only keep the variations that don't occur in NG c0 samples
  keepvar <- sapply(1:nrow(NGvariations),function(x) all(NGvariations[x,]==F))
  # apply keep to vcfR object
  vcf1 <- vcf[keepvar,]
  
  # remove variants where all samples have the same as the reference (artifact of removed problematic samples)
  filtempty <- data.frame(vcf1@gt)
  for (column in colnames(vcf1@gt[,-1])){
    filtempty[,column] <- startsWith(vcf1@gt[,column],"0:") | startsWith(vcf1@gt[,column],".:")
  }
  # get variants that can cause colistin resistance
  keepvar <- sapply(1:nrow(filtempty),function(x) any(filtempty[x,-1]==F))
  # apply keep to vcfR object
  vcf2 <- vcf1[keepvar,]
  
## potentional variants that cause colR ####
  gt <- extract.gt(vcf1, element = "GT",as.numeric = T)
  for (col in colnames(gt)){
    meta[which(meta$`new lable`==col),]$varCountFilt <- sum(gt[,col],na.rm = T)
  }
  
  # create dataframe for all possible colR variants in the strain
  if (strain == "KPN9884"){strainR <- "KPN9596-R"}
  if (strain == "KPN9749"){strainR <- "KPN9497-R"}
  
  sampleDF <- data.frame()
  for (sample in colnames(vcf2@gt[,-1])){
    # get vcf for the sample, select it's variations and get the annotation for the variations
    tempvcf <- vcf2[,c("FORMAT", sample)]
    pick <- startsWith(tempvcf@gt[,sample],"1:")
    if (length(which(pick))==0){next()}
    
    info <- extract.info(tempvcf, element = "ANN",as.numeric = F)
    
    # create a dataframe for all variations of the sample & add info
    variations <- data.frame("sample"=sample, "POS" = tempvcf@fix[pick,c("POS")], "REF"=tempvcf@fix[pick,c("REF")], "ALT"=tempvcf@fix[pick,c("ALT")] )
    variations$info <- info[pick]
    
    # bind the variationsDF to the sampleDF, creating one DF with all variants for the strain
    sampleDF <- rbind(sampleDF, variations)
    rownames(sampleDF) <- NULL
    
  }
  
  sampleDF <- merge(meta[,c("strain","replicate","Time point","CONC","new lable")], sampleDF, by.y="sample", by.x="new lable")
  
  # split sampleDF$info into multiple columns save into strainDF
  annotationCols <- c("ALT.info","annotation","putative_impact","gene_name","gene_ID","feature_type","feature_ID",
                      "Transcript_biotype","total","HGVS.c","HGVS.p","cDNA_position","CDS_position","Protein_position")
  
  strainDF <- data.frame(separate(sampleDF, info, into=annotationCols ,sep="\\|"))
  strainDF <- strainDF[!duplicated(strainDF),]
  strainDF$ALT.info <- NULL
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@testing@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  strainDF$gene_name
  ID <- "_00001"
  gff_line <- system(paste0("cat genome_data/toolOut/prokka/",strain,"/",strain,".gff | egrep ",ID),intern = T)
  
  product <- sub("^.*;product=", "", gff_line)
  inference <- sub(";.*","",sub("^.*;inference=", "", gff_line))
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@testing@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  
  

  
  
  # add dataframe of strain to the list of dataframes
  strainL[[strain]] <- strainDF
  
  
## create FASTA of hypothetical proteins #### 
  # get list of sequences for the hypothetical proteins
  hypo_prots <- strainDF[contains(hyponames[strain],vars = strainDF$gene_name),c("gene_name","gene_ID", "new.lable")]
  hypo_prots <- unique(hypo_prots$gene_ID)
  
  # read Nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
  sequences <- readDNAStringSet(paste0("genome_data/toolOut/prokka/",strain,"/",strain,".ffn"))
  
  # set file save name
  file_s <- "genome_data/toolOut/variantComparison_output/hypothetical_proteins.fa"
  
  for (protein in hypo_prots){
    if (grepl("-", protein)){ # two gene names
      write.table(paste0("##  ~~>",protein,"|",strain), file = file_s, append = T,row.names = F, quote = F,col.names = F) #write protein names to .fa file
      protein_u <- unlist(unique(strsplit(protein, "-")))
      for (i in protein_u){
        write.table(paste0(">",i,"|",strain), file = file_s, append = T, quote = F, row.names = F,col.names = F)
        prot_seq <- sequences[grep(i, names(sequences))]
        write.table(prot_seq, file = file_s, append = T, quote = F, row.names = F,col.names = F)
      }
      write.table(paste0("##  <~~",protein,"|",strain), file = file_s, append = T,row.names = F, quote = F,col.names = F)
    }else{ # one gene name
      write.table(paste0(">",protein,"|",strain), file = file_s, append = T,row.names = F, quote = F,col.names = F) #write protein names to .fa file
      prot_seq <- sequences[grep(protein, names(sequences))]
      write.table(prot_seq, file = file_s, append = T, quote = F, row.names = F,col.names = F)
    }
  }
}

# write variation dataframes to xlsx file
write_xlsx(strainL, "genome_data/toolOut/variantComparison_output/snpList.xlsx")

# get important columns and write to csv file
meta_save <- na.omit(meta[,c("strain", "experiment", "Time point", "CONC", "replicate", "New ID", "new lable","varCountRaw", "varCountFilt")])
write.csv(meta_save, file = "genome_data/toolOut/variantComparison_output/numberOfSNPs.csv")

quit()








#https://cran.rstudio.com/web/packages/g3viz/vignettes/introduction.html#31_example_1:_visualize_genetic_mutation_data_from_maf_file
  
  
  
  
  
  # create plots ###EXPERIMENTAL
  #vcfR
  gff <- read.delim(paste0("genome_data/toolOut/prokka/",strain,"/",strain,".gff"), header = F, comment.char = "#")
  gff.genes <- gff[gff[,3]=="gene",]
  
  dna <- ape::read.dna(gzfile("genome_data/toolOut/consensus/DTK_35_consensus.fa"), format = "fasta")
  
  mat.dna <- as.matrix(dna)
  
  chrom <- create.chromR(vcf=vcf1, ann = gff.genes, seq = mat.dna)
  gt <- extract.gt(vcf1, element = "GT",as.numeric = T)
  barplot(apply(gt, MARGIN=2, count, na.rm=TRUE),las = 3, main = strain)
  
  plot(chrom)
  
  chrom <- proc.chromR(chrom, verbose = T)
  chromoqc(chrom, dp.alpha = 22)
  
  
  
  # GenVisR
  library(GenVisR)
  
  waterfall(vcf1) 
  
  
  # variantannotation
  library(VariantAnnotation)
  
  
  
  
  
}


















##edit gewerkte uren
library(xlsx)
library(dplyr)

uren <- read_xlsx("gewerkte uren_n.xlsx") %>%
  mutate(
    datum = as.Date(as.numeric(datum), origin = "1899-12-30"), 
    begin = format(begin, "%H:%M"),   # “That’s what I do: I drink and I know things.”
    eind = format(eind, "%H:%M"),
    totaal = format(totaal, "%H:%M"),
    `totaal week` = format(`totaal week`, "%H:%M")
)


as.Date(uren$begin, origin = "1/1/1990 00:00")
uren[68,c("begin","eind","totaal")] <- c("11:00","14:00","3:00")
uren[67,c("begin","eind","totaal")] <- c("10:30","17:30","7:00")

write_xlsx(uren, "gewerkte uren_n.xlsx")


















  
  