# plotting script is inspired by Linda Aulin

## Clean up workspace
rm(list=ls())

##  LOAD PACKAGES & USER AGRUMENTS  ###################################################################################################################
pkgs <- c("dplyr", "readxl","tidyr","ggplot2","ggrepel")
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
outputdir <- args[1]
inputfile <- args[2]



##  PLOT DATA  ###################################################################################################################
plot_data = function(dat, strain){
  new_dat <- dat %>%
    separate(col = new.lable,into = c("Experiment", "ID"), sep = "_") %>% 
    mutate(TIME = gsub("p2", 48, Time.point)) %>%
    mutate(TIME = as.numeric(gsub("t", "", TIME)))
  
  # Filter out any duplicated rows based on gene name, strain, TIME, replicate, Experiment for plotting
  final_dat <- new_dat %>% 
    distinct(gene_name, strain, TIME, replicate, Experiment, .keep_all = T)
  
  ggplot(final_dat, aes(y = Experiment, x= gene_name))+
    geom_point(aes( group = replicate, shape = replicate), size = 5, alpha = 0.5, position = position_dodge(width = 0.5))+
    facet_grid(~TIME)+
    scale_alpha_continuous(range = c(0.5,1)  )+
    labs(y = "Experiment",
         x = "Gene name")+
    coord_flip()+
    theme_bw()+
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle(strain)
}

# read excel file and plot with plot_data for each sheet, save to plots
plots <- list()
for (i in 1:length( excel_sheets( paste0(outputdir,inputfile) ) )){
  
  # get strain name (sheetname)
  strain <- excel_sheets("genome_data/toolOut/results/2020-02-25_snpList.xlsx")[i]
  
  # create plot with plot_data function
  plots[[i]] <- plot_data(read_excel(paste0(outputdir,inputfile), sheet= i), strain)
}


# write plots to .pdf
pdf(file = paste0(outputdir,sub(".xlsx","",inputfile),".pdf"), width = 10, height = 10)
for (figure in plots){
  plot(figure)
}
dev.off()

quit()