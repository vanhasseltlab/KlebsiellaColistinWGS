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
clear
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
# guideline for chosing the right kmer size, please try multiple sizes and visualize (e.g. bandage)
# https://github.com/rrwick/Bandage/wiki/Effect-of-kmer-size
  [ ! -s ${outdir}/shovill/${strain}/contigs.fa ] &&\
  shovill --R1 $R1 --R2 $R2 --outdir ${outdir}/shovill/${strain}/ \
  --kmers 101 --mincov 10 --minlen 500 --force\
  --assembler spades --opts --careful --tmpdir /tmp

##~~~~ prokka: annotate consensus.fa reference  ~~~~##
  [ ! -d ${outdir}/prokka/${strain} ] &&\
  prokka --force --addgenes --prefix $strain --outdir ${outdir}/prokka/${strain} ${outdir}/shovill/${strain}/contigs.fa

##~~~~ snpEff: build database for baseline contigs.fa references  ~~~~##
if [ ! -f /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil/genes.gbk ];then
  # create needed directory
  mkdir -p /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil
  
  # copy genbank file and gff file for creating the database
  cp ${outdir}/prokka/${strain}/$strain.gbk /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil
  mv /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil/$strain.gbk /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil/genes.gbk
  
  cp ${outdir}/prokka/${strain}/$strain.gff /home/linuxbrew/.linuxbrew/share/snpeff/data/${strain}_shovil
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






###################### testing multiple assembly settings
# add to PATH for Prokka
PATH=${PATH}:/home/linuxbrew/.linuxbrew/bin/

Reads="genome_data/raw_sequences/"
label="DTK_01"
label="DTK_18"

R1=$(ls ${Reads}*R1* | egrep ${label})
R2=$(echo ${R1} | sed 's/R1.filt.fastq.gz$/R2.filt.fastq.gz/g')

mkdir -p test_assemblies

# kmer comparisons
for kmer in 77 99 101 105;do
  echo $kmer
    shovill --R1 $R1 --R2 $R2 --outdir test_assemblies/shovil_spades_kmer${kmer}/ \
  --kmers ${kmer} --mincov 10 --minlen 200 --force\
  --assembler spades

  shovill --R1 $R1 --R2 $R2 --outdir test_assemblies/shovil_spades-careful_kmer${kmer}/ \
  --kmers ${kmer} --mincov 10 --minlen 200 --force\
  --assembler spades --opts --careful
done
  
tools/quast-5.0.2/quast.py -o test_assemblies/results_kmer_comparison\
test_assemblies/shovil_spades_kmer77/contigs.fa test_assemblies/shovil_spades-careful_kmer77/contigs.fa\
test_assemblies/shovil_spades_kmer99/contigs.fa test_assemblies/shovil_spades-careful_kmer99/contigs.fa\
test_assemblies/shovil_spades_kmer101/contigs.fa test_assemblies/shovil_spades-careful_kmer101/contigs.fa\
test_assemblies/shovil_spades_kmer105/contigs.fa test_assemblies/shovil_spades-careful_kmer105/contigs.fa
  

# tool comparisons
for assembler in spades skesa velvet megahit;do
  echo $assembler
  shovill --R1 $R1 --R2 $R2 --outdir test_assemblies/shovil_{assembler}_kmer101/ \
  --kmers 101 --mincov 10 --minlen 200 --force\
  --assembler $assembler
done
  
tools/quast-5.0.2/quast.py -o test_assemblies/results_tool_comparison\
test_assemblies/shovil_spades_kmer101/contigs.fa\
test_assemblies/shovil_skesa_kmer101/contigs.fa\
test_assemblies/shovil_velvet_kmer101/contigs.fa\
test_assemblies/shovil_megahit_kmer101/contigs.fa
###################### testing multiple settings
  
