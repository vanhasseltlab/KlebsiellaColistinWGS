#!/bash/bin
#
#===========================================================================================================#
#                           Introduction
# 
# This script does De Novo assembly with shovill from the linuxbrew package.
# Then it annotates the contigs.fa file with Prokka.
# From which snpEff databases are build for the annotation of variants after variant calling.
#
usage(){
echo "
deNovoAssembler version 0.0.1 part of the runNGSdownstreamPipeline (build 2020-02-18), by Yob Haakman <yob1997@live.nl>.

Usage:  bash deNovoAssembler.sh [STRAIN] [PATH_TO_R1.fastq.gz] [PATH_TO_R2.fastq.gz] [OUTPUT_DIRECTORY]

Required parameters:
  [STRAIN]                :Strain name, will be used in annotation and for building a snpEff database
  [PATH_TO_R1.fastq.gz]   :Path to forward read of sample (fastq.gz format)
  [PATH_TO_R2.fastq.gz]   :Path to reverse read of sample (fastq.gz format)
  [OUTPUT_DIRECTORY]      :Folder to output results to [stdout]

Input parameters have to be in the order displayed as above.
"
}

if [ -z "${1}" ] || [ "${1}" == "-h" ] || [ "${1}" == "-help" ];then
  usage
  exit
fi
#===========================================================================================================#



#===========================================================================================================#
#                           Set variables

strain=$1; R1=$2; R2=$3; outdir=$4
echo $strain $R1 $R2 $outdir
#===========================================================================================================#



#===========================================================================================================#
#                           De novo assembly, reordering and cleaning

##~~~~ Shovil: de novo assembly and refinementwith best tested settings ~~~~##    
  [ ! -s ${outdir}/shovill/${strain}/contigs.fa ] &&\
  shovill --R1 $R1 --R2 $R2 --outdir ${outdir}/shovill/${strain}/ \
  --kmers 101 --mincov 10 --minlen 500 --force\
  --assembler spades --opts --careful
  
##~~~~ prokka: annotate consensus.fa reference  ~~~~##
  [ ! -d ${outdir}/prokka/${strain} ] &&\
  prokka --force --addgenes --prefix $strain --outdir ${outdir}/prokka/${strain} ${outdir}/shovill/${strain}/contigs.fa

##~~~~ snpEff: build database for baseline contigs.fa references  ~~~~##
if [ ! -f /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil/genes.gbk ];then
  # create needed directory
  mkdir -p /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil
  
  # copy genbank file and gff file for creating the database
  cp genome_data/toolOut/prokka/${strain}/$strain.gbk /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil
  mv /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil/$strain.gbk /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil/genes.gbk
  
  cp genome_data/toolOut/prokka/${strain}/$strain.gff /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil
  mv /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil/$strain.gff /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil/regulation.gff
  
  # add database to config file, if it isn't in there yet
  if ! grep -Fxq "${strain}_shovil.genome : ${strain}_shovil" /home/linuxbrew/.linuxbrew/share/snpeff/snpEff.config;then
    echo "${strain}_shovil.genome : ${strain}_shovil" >> /home/linuxbrew/.linuxbrew/share/snpeff/snpEff.config
  fi
  
  # build snpEff database
  snpEff build -genbank -v ${strain}_shovil
fi  
exit
#                           End of script
#===========================================================================================================#


















# add to PATH for Prokka
PATH=${PATH}:/home/linuxbrew/.linuxbrew/bin/

# Rename the numeric variables to their original names, for better script readabiliy.
Reads="genome_data/raw_sequences/"
assembly="genome_data/toolOut/SPAdes_assembly/"
redundans="genome_data/toolOut/redundans/"
consensus="genome_data/toolOut/consensus/"

# and make directories if they dont exist
mkdir -p ${assembly} ${consensus} ${redundans}

declare -A strains=( ["DTK_01"]="KPN9884" ["DTK_18"]="KPN9749" )
#===========================================================================================================#



#===========================================================================================================#
#                           De novo assembly, reordering and cleaning

for label in ${!strains[*]};do
  R1=$(ls ${Reads}*R1* | egrep ${label})
  R2=$(echo ${R1} | sed 's/R1.filt.fastq.gz$/R2.filt.fastq.gz/g')
  strain=${strains[$label]}
  
  echo $label ${strains[$label]} $R1 $R2
  echo ""

##~~~~ SPAdes: de novo assembly with best settings tested ~~~~##    
  [ ! -s ${assembly}${strain}_shovil/contigs.fa ] &&\
  shovill --R1 $R1 --R2 $R2 --outdir ${assembly}${strain}_shovil/ \
  --kmers 101 --mincov 10 --minlen 200 --force\
  --assembler spades --opts --careful
  
  
##~~~~ prokka: annotate consensus.fa reference  ~~~~##
  [ ! -d genome_data/toolOut/prokka/${strain} ] &&\
  prokka --force --addgenes --prefix $strain --outdir genome_data/toolOut/prokka/${strain} ${assembly}${strain}_shovil/contigs.fa

##~~~~ snpEff: build database for baseline contigs.fa references  ~~~~##
if [ ! -f /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}/genes.gbk ];then
  # create needed directory
  mkdir -p /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil
  
  # copy genbank file and gff file for creating the database
  cp genome_data/toolOut/prokka/${strain}/$strain.gbk /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil
  mv /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil/$strain.gbk /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil/genes.gbk
  
  cp genome_data/toolOut/prokka/${strain}/$strain.gff /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil
  mv /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil/$strain.gff /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil/regulation.gff
  
  # add database to config file, if it isn't in there yet
  if ! grep -Fxq "${strain}_shovil.genome : ${strain}_shovil" /home/linuxbrew/.linuxbrew/share/snpeff/snpEff.config;then
    echo "${strain}_shovil.genome : ${strain}" >> /home/linuxbrew/.linuxbrew/share/snpeff/snpEff.config
  fi
  
  # build snpEff database
  snpEff build -genbank -v ${strain}_shovil
fi  
done
#===========================================================================================================#
exit






declare -A references=( ["DTK_01"]="CP034281.1" ["DTK_18"]="CP031795.1" )

# improve the de novo assembly by ordering the contigs.fa on the ncbi template
# and create scaffolds and fill gaps with redundands. The scaffolds.filled.fa is the final assembly
  [ ! -d ${redundans}${label} ] &&\
  tools/redundans/redundans.py \
  -i $R1 $R2 \
  -r ${assembly}${label}_NCBItemplate.fa \
  -f ${assembly}${label}/contigs.fasta \
  -o ${redundans}${label}/ -t 16 --noscaffolding
  
  
  [ ! -s ${consensus}${baseline}_NCBItemplate.fa ] &&\
  esearch -db nucleotide -query ${references[$label]} | efetch -db nucleotide -format fasta > ${assembly}${label}_NCBItemplate.fa
  



###################### testing multiple settings
  [ ! -s ${assembly}${label}_deNovo/contigs.fasta ] &&\
  tools/SPAdes-3.13.0-Linux/bin/spades.py -1 ${R1} -2 ${R2} -o ${assembly}${label} -t 32 -k 77,99,101,105 &&\
  tools/SPAdes-3.13.0-Linux/bin/spades.py -1 ${R1} -2 ${R2} -o ${assembly}${label} -t 32 -k 77,99,101,105 --careful

  tools/quast-5.0.2/quast.py \
  genome_data/toolOut/shovil_kmer77/contigs.fa genome_data/toolOut/shovil_kmer77_careful/contigs.fa \
  genome_data/toolOut/shovil_kmer99/contigs.fa genome_data/toolOut/shovil_kmer99_careful/contigs.fa \
  genome_data/toolOut/shovil_kmer101/contigs.fa genome_data/toolOut/shovil_kmer101_careful/contigs.fa \
  genome_data/toolOut/shovil_kmer105/contigs.fa genome_data/toolOut/shovil_kmer105_careful/contigs.fa
  
  shovill --R1 $R1 --R2 $R2 --outdir genome_data/toolOut/shovil_spades-careful_kmer101/ \
  --kmers 101 --mincov 10 --minlen 200 --force\
  --assembler spades --opts --careful

  shovill --R1 $R1 --R2 $R2 --outdir genome_data/toolOut/shovil_skesa_kmer101/ \
  --kmers 101 --mincov 10 --minlen 200 --force\
  --assembler skesa 
  
  shovill --R1 $R1 --R2 $R2 --outdir genome_data/toolOut/shovil_velvet_kmer101/ \
  --kmers 101 --mincov 10 --minlen 200 --force\
  --assembler velvet
  
  shovill --R1 $R1 --R2 $R2 --outdir genome_data/toolOut/shovil_megahit_kmer101/ \
  --kmers 101 --mincov 10 --minlen 200 --force\
  --assembler megahit 
  
  tools/quast-5.0.2/quast.py -o quast_results/results_tool_comparison\
  genome_data/toolOut/shovil_spades-careful_kmer101/contigs.fa \
  genome_data/toolOut/shovil_skesa_kmer101/contigs.fa \
  genome_data/toolOut/shovil_velvet_kmer101/contigs.fa \
  genome_data/toolOut/shovil_megahit_kmer101/contigs.fa
  
  
  ###################### testing multiple settings
  
