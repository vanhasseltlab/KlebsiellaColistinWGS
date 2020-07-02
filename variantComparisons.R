## Clean up workspace
rm(list=ls())



##  LOAD PACKAGES & USER AGRUMENTS  ###################################################################################################################
# packages
pkgs <- c("vcfR", "readxl","writexl","tidyr","stringr")
for (pkg in pkgs){
  if (!pkg %in% installed.packages()) {
    install.packages(pkg, dependencies = TRUE)
  }
}
lapply(pkgs, library, character.only = TRUE)
remove (pkgs,pkg)

## load user arguments
args <- commandArgs(T)
# output directory
outdir <- args[1]
dir.create(file.path(outdir, "results"), showWarnings = FALSE)

# path to meta files
meta_file <- args[2]

outdir <- "genome_data/toolOut/"
meta_file <- "Sample_overview_6-23-2020.xlsx"


##  READ DATA  #########################################################################################################################
# meta data of samples
meta <- read_excel(meta_file)
meta <- as.data.frame(meta)



##  RUN ANALYSIS  #########################################################################################################################
# create lists and dataframes to save data to
strainL <- list()
strainDF <- data.frame()
strainL_noNG <- list()
strainDF_noNG <- data.frame()

# loop through unique strains in meta
for (strain in unique(meta$strain)){
  
  # read vcf file with vcfR
  vcf <- read.vcfR(paste0(outdir,"/GATK_pipe_output/",strain,"/annotated.vcf.gz"))
  #vcf <- read.vcfR(paste0("pgap_run/KPN9749_annotated.vcf.gz"))
  
  ## WITHOUT NG CLEANING ####
  sampleDF_noNG <- data.frame()
  for (sample in colnames(vcf@gt[,-1])){
    # get vcf for the sample, select it's variations and get the annotation for the variations
    tempvcf <- vcf[,c("FORMAT", sample)]
    pick <- startsWith(tempvcf@gt[,sample],"1:")
    if (length(which(pick))==0){next()}
    
    info <- extract.info(tempvcf, element = "ANN",as.numeric = F)
    
    # create a dataframe for all variations of the sample & add info
    variations <- data.frame("sample"=sample,"contig"=tempvcf@fix[pick,c("CHROM")],"POS" = tempvcf@fix[pick,c("POS")], "REF"=tempvcf@fix[pick,c("REF")], "ALT"=tempvcf@fix[pick,c("ALT")] )
    variations$info <- info[pick]
    
    # bind the variationsDF to the sampleDF_noNG, creating one DF with all variants for the strain
    sampleDF_noNG <- rbind(sampleDF_noNG, variations)
    rownames(sampleDF_noNG) <- NULL
  }
  
  ## WITH NG CLEANING  ####
  # get all NG c0 samples and put their genotypes in a dataframe.
  NGsamples <- meta[which(meta$replicate=="NG" & meta$strain==strain & meta$CONC=="c0"),]$`new lable`
  NGvariations <- data.frame(vcf@gt[,c(colnames(vcf@gt)%in%NGsamples)])
  # loop though the columns of the NG dataframe, check if the NG samples have the variation (TRUE/FALSE)
  for (column in NGsamples){
    print (column)
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
  
  ## potentional variants that cause colR
  sampleDF <- data.frame()
  for (sample in colnames(vcf2@gt[,-1])){
    # get vcf for the sample, select it's variations and get the annotation for the variations
    tempvcf <- vcf2[,c("FORMAT", sample)]
    pick <- startsWith(tempvcf@gt[,sample],"1:")
    if (length(which(pick))==0){next()}
    
    info <- extract.info(tempvcf, element = "ANN",as.numeric = F)
    
    # create a dataframe for all variations of the sample & add info
    variations <- data.frame("sample"=sample,"contig"=tempvcf@fix[pick,c("CHROM")], "POS" = tempvcf@fix[pick,c("POS")], "REF"=tempvcf@fix[pick,c("REF")], "ALT"=tempvcf@fix[pick,c("ALT")] )
    variations$info <- info[pick]
    
    # bind the variationsDF to the sampleDF, creating one DF with all variants for the strain
    sampleDF <- rbind(sampleDF, variations)
    rownames(sampleDF) <- NULL
  }
  
  
  ## ADDITIONAL ANNOTATION  ######################################################################################################################### 
  sampleDF <- merge(meta[,c("strain","replicate","Time point","CONC","new lable")], sampleDF, by.y="sample", by.x="new lable")
  sampleDF_noNG <- merge(meta[,c("strain","replicate","Time point","CONC","new lable")], sampleDF_noNG, by.y="sample", by.x="new lable")
  
  # split sampleDF$info into multiple columns save into strainDF
  annotationCols <- c("ALT.info","annotation","putative_impact","gene_name","gene_ID","feature_type","feature_ID",
                      "Transcript_biotype","total","HGVS.c","HGVS.p","cDNA_position","CDS_position","Protein_position")
  
  strainDF <- data.frame(separate(sampleDF, info, into=annotationCols ,sep="\\|"))
  strainDF <- strainDF[!duplicated(strainDF),]
  strainDF$ALT.info <- NULL
  
  strainDF_noNG <- data.frame(separate(sampleDF_noNG, info, into=annotationCols ,sep="\\|"))
  strainDF_noNG <- strainDF_noNG[!duplicated(strainDF_noNG),]
  strainDF_noNG$ALT.info <- NULL
  
  # add contig length to dataframes
  strainDF$contig_length <- NA
  for (cont in unique(strainDF$contig)){
    strainDF[which(strainDF$contig==cont),]$contig_length <- gsub(".*length=(.+),.*", "\\1", grep(cont,vcf@meta,value=T))
  }
  
  strainDF_noNG$contig_length <- NA
  for (cont in unique(strainDF_noNG$contig)){
    strainDF_noNG[which(strainDF_noNG$contig==cont),]$contig_length <- gsub(".*length=(.+),.*", "\\1", grep(cont,vcf@meta,value=T))
  }
  
  
  # add additional annotation from gff files 
  strainDF$product <- NA
  strainDF$inference <- NA
  for (ID in unique(grep("^pgaptmp_[0-9]+.*$",strainDF$gene_ID, value=T))){
    
    if (grepl("-",ID)){
      IDn <- sub("[A-Z]+_","",unlist(strsplit(ID, "-")))
      #gff_line <- system(paste0("cat ",outdir,"/prokka/",strain,"/",strain,".gff | egrep '",IDn[1],"|",IDn[2],"'"),intern = T)
      gff_line <- system(paste0("cat pgap_run/output_",strain,"/annot.gff | egrep '",IDn[1],"|",IDn[2],"'"),intern = T)
      
      #product <- paste(sub("^.*;product=", "", gff_line[grep("product",gff_line)]),collapse = " ;-; ")
      product <- paste(gsub(".*product=(.+);p.*", "\\1", gff_line[grep("product",gff_line)]),collapse = " ;-; ")
      inference <- paste(gsub(".*inference=(.+);l.*", "\\1", gff_line[grep("inference",gff_line)]),collapse = " ;-; ")
    }else{
      IDn <- sub("[A-Z]+_","",ID)
      #gff_line <- system(paste0("cat ",outdir,"/prokka/",strain,"/",strain,".gff | egrep ",IDn),intern = T)
      gff_line <- system(paste0("cat pgap_run/output_",strain,"/annot.gff | egrep ",IDn),intern = T)
      
      
      #product <- sub("^.*;product=", "", gff_line[grep("product",gff_line)])
      product <- gsub(".*product=(.+);p.*", "\\1", gff_line)
      inference <- gsub(".*inference=(.+);l.*", "\\1", gff_line)
    }
    strainDF[which(strainDF$gene_ID==ID),]$product <- product[2]
    strainDF[which(strainDF$gene_ID==ID),]$inference <- inference[2]
  }
  
  strainDF_noNG$product <- NA
  strainDF_noNG$inference <- NA
  for (ID in unique(grep("^pgaptmp_[0-9]+.*$",strainDF_noNG$gene_ID, value=T))){
    
    if (grepl("-",ID)){
      IDn <- sub("[A-Z]+_","",unlist(strsplit(ID, "-")))
      #gff_line <- system(paste0("cat ",outdir,"/prokka/",strain,"/",strain,".gff | egrep '",IDn[1],"|",IDn[2],"'"),intern = T)
      gff_line <- system(paste0("cat pgap_run/output_",strain,"/annot.gff | egrep '",IDn[1],"|",IDn[2],"'"),intern = T)
      #product <- paste(sub("^.*;product=", "", gff_line[grep("product",gff_line)]),collapse = " ;-; ")
      product <- paste(gsub(".*product=(.+);p.*", "\\1", gff_line[grep("product",gff_line)]),collapse = " ;-; ")
      inference <- paste(gsub(".*inference=(.+);l.*", "\\1", gff_line[grep("inference",gff_line)]),collapse = " ;-; ")
    }else{
      IDn <- sub("[A-Z]+_","",ID)
      #gff_line <- system(paste0("cat ",outdir,"/prokka/",strain,"/",strain,".gff | egrep ",IDn),intern = T)
      gff_line <- system(paste0("cat pgap_run/output_",strain,"/annot.gff | egrep ",IDn),intern = T)
      product <- gsub(".*product=(.+);p.*", "\\1", gff_line)
      inference <- gsub(".*inference=(.+);l.*", "\\1", gff_line)
    }
    strainDF_noNG[which(strainDF_noNG$gene_ID==ID),]$product <- product[2]
    strainDF_noNG[which(strainDF_noNG$gene_ID==ID),]$inference <- inference[2]
  }
  
  
  # some snps have alt=* these don't have the correct gene_ID, gene_name, gene_feature tag, 
  #this is fixed by reading the gff file and comparing the position of the snp to the genes 
  contigs <- unique(strainDF[which(strainDF$ALT %in% "*"),]$contig)
  for (contig in contigs){
    skip_to_next <- FALSE
    tryCatch({
      gff_lines <- system(paste0("cat pgap_run/output_",strain,"/annot.gff | egrep ",contig,' | grep -e "CDS"'),intern = T)
      
      gff_lines <- data.frame(do.call('rbind', strsplit(as.character(gff_lines),'\t',fixed=TRUE)))
      
      gff_lines[,c("X4","X5")] <- sapply(gff_lines[,c("X4","X5")],function(x) as.numeric(as.character(x)))
      snp_positions <- as.numeric(as.character(strainDF[which(strainDF$ALT %in% "*"),]$POS))
      
      for (i in 1:length(gff_lines$X4)){
        test <- snp_positions[ snp_positions >= gff_lines$X4[i] & snp_positions <= gff_lines$X5[i]]
        if (length(test)!=0){
          gff_line <- gff_lines[i,]$X9
          strainDF[which(strainDF$POS %in% snp_positions),]$product <- gsub(".*product=(.+);p.*", "\\1", gff_line)
          strainDF[which(strainDF$POS %in% snp_positions),]$inference <- gsub(".*inference=(.+);l.*", "\\1", gff_line)
          strainDF[which(strainDF$POS %in% snp_positions),c("gene_name","gene_ID","feature_ID")] <- gsub(".*ID=(.+);Parent.*", "\\1", gff_line)
          
        }
      }
    }, error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next) { next } 
  }
  
  # contigs <- unique(strainDF_noNG[which(strainDF_noNG$ALT %in% "*"),]$contig)
  # for (contig in contigs){
  #   skip_to_next <- FALSE
  #   tryCatch({
  #     #gff_lines <- system(paste0("cat ",outdir,"/prokka/",strain,"/",strain,".gff | egrep ",contig,' | grep -e "CDS" -e "Aragorn"'),intern = T)
  #     gff_lines <- system(paste0("cat pgap_run/output_",strain,"/annot.gff | egrep ",contig,' | grep -e "CDS"'),intern = T)
  #   }, error = function(e) { skip_to_next <<- TRUE})
  #   if(skip_to_next) { next }
  #     
  #     gff_lines <- data.frame(do.call('rbind', strsplit(as.character(gff_lines),'\t',fixed=TRUE)))
  # 
  #     
  #     gff_lines[,c("X4","X5")] <- sapply(gff_lines[,c("X4","X5")],function(x) as.numeric(as.character(x)))
  #     snp_positions <- as.numeric(as.character(strainDF_noNG[which(strainDF_noNG$ALT %in% "*"),]$POS))
  #     
  #     for (i in 1:length(gff_lines$X4)){
  #       test <- snp_positions[ snp_positions >= gff_lines$X4[i] & snp_positions <= gff_lines$X5[i]]
  #       if (length(test)!=0){
  #         gff_line <- gff_lines[i,]$X9
  #         strainDF_noNG[which(strainDF_noNG$POS %in% snp_positions),]$product <- sub("^.*;product=", "", gff_line)
  #         strainDF_noNG[which(strainDF_noNG$POS %in% snp_positions),c("gene_name","gene_ID","feature_ID")] <- gsub(".*ID=(.+);Parent.*", "\\1", gff_line)
  #       }
  #     }
  # }
  
  # reorder dataframe, to a better readable order
  col_order <- c("strain", "Time.point", "replicate","CONC", "new.lable", "contig", "contig_length", "POS", "REF", "ALT",  "gene_name","product",
                 "inference", "HGVS.c", "HGVS.p", "cDNA_position", "CDS_position", "Protein_position", "annotation", "putative_impact","gene_ID", "feature_type", "feature_ID","Transcript_biotype", "total")
  strainDF <- strainDF[, col_order]
  strainDF_noNG <- strainDF_noNG[, col_order]
  
  # make gene_names of ORFs prettier
  for (ID in unique(grep("^.*pgaptmp_[0-9]+.*$",strainDF$gene_ID, value=T))){
    if (grepl("[0-9]-",ID)){
      IDn <- sub("[cds\\-]*pgaptmp_","",unlist(strsplit(ID, "-")))
      new_name <- paste0("ORF_",IDn[1],"-ORF_",IDn[2])
      try(strainDF[which(strainDF$gene_name==ID),]$gene_name <- new_name,silent=T)
      try(strainDF[which(strainDF$gene_ID==ID),]$gene_ID <- new_name,silent=T)
      try(strainDF[which(strainDF$feature_ID==ID),]$feature_ID <- new_name,silent=T)
    }else{
      IDn <- sub("[cds\\-]*pgaptmp_","",ID)
      new_name <- paste0("ORF_",IDn)
      try(strainDF[which(strainDF$gene_name==ID),]$gene_name <- new_name,silent=T)
      try(strainDF[which(strainDF$gene_ID==ID),]$gene_ID <- new_name,silent=T)
      try(strainDF[which(strainDF$feature_ID==ID),]$feature_ID <- new_name,silent=T)
    }
  }
  
  # make gene_names of ORFs prettier
  for (ID in unique(grep("^pgaptmp_[0-9]+.*$",strainDF_noNG$gene_ID, value=T))){
    if (grepl("[0-9]-",ID)){
      IDn <- sub("[cds\\-]*pgaptmp_","",unlist(strsplit(ID, "-")))
      new_name <- paste0("ORF_",IDn[1],"-ORF_",IDn[2])
      try(strainDF_noNG[which(strainDF_noNG$gene_name==ID),]$gene_name <- new_name,silent=T)
      try(strainDF_noNG[which(strainDF_noNG$gene_ID==ID),]$gene_ID <- new_name,silent=T)
      try(strainDF_noNG[which(strainDF_noNG$feature_ID==ID),]$feature_ID <- new_name,silent=T)
    }else{
      IDn <- sub("[cds\\-]*pgaptmp_","",ID)
      new_name <- paste0("ORF_",IDn)
      try(strainDF_noNG[which(strainDF_noNG$gene_name==ID),]$gene_name <- new_name,silent=T)
      try(strainDF_noNG[which(strainDF_noNG$gene_ID==ID),]$gene_ID <- new_name,silent=T)
      try(strainDF_noNG[which(strainDF_noNG$feature_ID==ID),]$feature_ID <- new_name,silent=T)
    }
  }
  
  # add dataframe of strain to the list of dataframes
  strainL[[strain]] <- strainDF
  strainL_noNG[[strain]] <- strainDF_noNG
  
  ## create FASTA of hypothetical proteins (turned off) #### 
  if (FALSE){
    # get list of sequences for the hypothetical proteins
    hypo_prots <- unique(strainDF[contains("hypothetical protein",vars = strainDF$product),c("gene_ID")])
    
    # read Nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
    sequences <- readDNAStringSet(paste0(outdir,"/prokka/",strain,"/",strain,".ffn"))
    
    # set file save name
    file_s <- paste0(outdir,"/hypothetical_proteins.fa")
    
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
}

# write output xlsx files
write_xlsx(strainL, paste0(outdir,"/results/",Sys.Date(),"_NG-filteration_snpList.xlsx"))
write_xlsx(strainL_noNG, paste0(outdir,"/results/",Sys.Date(),"_snpList.xlsx"))

#write_xlsx(strainL, paste0("pgap_run/",Sys.Date(),"_KPN9749_snpList.xlsx"))


quit()