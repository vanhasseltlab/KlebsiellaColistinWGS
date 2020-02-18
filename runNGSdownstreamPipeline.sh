#!/bash/bin
#===========================================================================================================#
#                           Introduction

usage(){
echo "
runNGSdownstramPipeline version 0.0.1 (build 2020-02-17), by Yob Haakman <yob1997@live.nl>.
This script does variation calling on filtered Next Generation Sequencing files (fastq.gz format)
Creates a draft genome from the inoculum sample reads, annotates the draft genome and create databases for it.
Maps evolved isogenic mutant sample reads to the draft genome.
Does variation calling based on the draft genome.
Output can be used in downstream analysis with R.

Usage:  bash runNGSdownstreamPipeline.sh -r [PATH_TO_R1_AND_R2_READS] -m [PATH_TO_META.xlsx]

Required parameters:
  -r , -reads <path>      :Path to reads.
                          Input must be a folder containing forward (R1) and reverse (R2) reads for each sample in compressed fastq format (fastq.gz).
  -m , -meta <path>       :Path to meta xlsx file. (default: Sample_overview_seq.xlsx)
                          (Tab-delimited xlsx file with columns in the following order.
                          1:'strain'  2:'experiment' 3:'typ of sample'  4:'Time point'	5:'CONC'	6:'replicate'	
                          7:'experimentalist'	8:'date'	9:'New ID'	10:'new lable'	11:'Box'	12:'row'	13:'col'
                          column 10 has to be unique for every sample.)

Optional parameters:
-o , -output <path>      :This will set the output folder to which all files are writen to. (default: NGS_pipe_output)
"
}

if [ -z "${1}" ] || [ "${1}" == "-h" ] || [ "${1}" == "-help" ];then
  usage
  exit
fi
#===========================================================================================================#



#===========================================================================================================#
#                           Check dependencies

# add to linuxbrew bin/ to PATH
PATH=${PATH}:/home/linuxbrew/.linuxbrew/bin/

for tool in "shovill" "prokka" "snpEff" "SnpSift" "vt" "tabix" "bgzip" "bwa mem" "xlsx2csv";do
  command -v $tool >/dev/null 2>1 || { echo >&2 "I require $tool but it's not installed.  Aborting."; exit 1; }
done

if [ ! -s tools/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar ];then
  echo "tools/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar not found.  Aborting."; exit 1;
fi

if [ ! -s tools/picard.jar ];then
  echo "tools/picard.jar not found.  Aborting."; exit 1;
fi
#===========================================================================================================#



#===========================================================================================================#
#                           Get input info

while test $# -gt 0; do # loop over input arguments
  case "$1" in
  -o|-output)
    shift
    outdir=$1
    shift
    ;;
  -r|-reads)
    shift
    reads=$1
    shift
    ;;
  -m|-meta)
    shift
    meta=$1
    shift
    ;; 
  *)
    echo "$1 is not a recognized option!"
  exit 1;
  ;;
  esac
done  

# check if reads option is given, directory exists and is non-empty
if [ -z "${reads}" ];then 
  echo "Please set option -reads [PATH_TO_READS].   Aborting."
  exit
elif [ ! -d "${reads}" ];then
  echo "Reads directory doesn't exist.   Aborting."
  exit
elif [ -z "$(ls -A ${reads})" ];then
  echo "Reads directory is empty.   Aborting."
  exit
fi

outdir=${outdir:-NGS_pipe_output/}
mkdir -p $outdir

# read meta file, fill array with variables
meta=${meta:-'Sample_overview_seq.xlsx'}
List_of_strains=$(xlsx2csv -i --skipemptycolumns -d "\t" "${meta}" | tail -n +2  | awk '{print $1}' | sort -u)

declare -A strains
for strain in $List_of_strains;do
  label=$(xlsx2csv -i -d "\t" "${meta}" | egrep -E "${strain}.*t0.*inoc" | sort -u | awk '{print $10}')
  strains[${label}]=${strain}
done

# check if array is not empty, and every strain has a key
if [ ${#strains[@]} -eq 0 ];then
  echo "Something is wrong with the meta xlsx file.   Aborting"; exit 1;
fi
for key in ${!strains[@]};do
  if [ -z ${strains[$key]} ];then
    echo "Something is wrong with the meta xlsx file.   Aborting"; exit 1;
  fi
done

echo "
Starting analysis NGS!

Output directory : ${outdir}
Path to reads : ${reads}
Path to meta file: ${meta}
"

#===========================================================================================================#



#===========================================================================================================#
#                           Create draft genomes of inoculum

for baseline in "${!strains[@]}";do
  (
  if [ ! -f ${outdir}/shovill/${strains[$baseline]}/contigs.fa ];then
    R1=$(ls ${reads}*R1* | egrep ${baseline})
    R2=$(echo ${R1} | sed 's/R1.filt.fastq.gz$/R2.filt.fastq.gz/g')
    bash DeNovoVarCalling/deNovoAssembler.sh ${strains[$baseline]} $R1 $R2 $outdir
    echo ""
  fi
  )&
done
wait


clear
echo "Baseline references created!"
echo "Next steps can take a long time to run."

while true;do
 read -p "Continue with pipeline? [y/n]" input
 case $input in
     [yY][eE][sS]|[yY])
 break;;
     [nN][oO]|[nN])
 exit;;
     *)
 echo "Invalid input...";;
 esac
done
#===========================================================================================================#



#===========================================================================================================#
#                           Mapping to contigs.fa from the deNovoassembler.sh

for baseline in "${!strains[@]}";do
  echo $baseline

  R1=$(ls ${reads}*R1* | egrep ${baseline})                             #Path to forward reads
  R2=$(echo ${R1} | sed 's/R1.filt.fastq.gz$/R2.filt.fastq.gz/g')       #Path to reverse reads
  ref="${outdir}/shovill/${strains[$baseline]}/contigs.fa"              #Path to reference consensus genome
  
##~~~~ BWA MEM/GATK/samtools: index reference genome ~~~~##
  [ ! -s $ref.ann ] && bwa index $ref
  [ ! -s ${outdir}/shovill/${strains[$baseline]}/contigs.dict ] && java -jar tools/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar CreateSequenceDictionary -R $ref
  [ ! -s $ref.fai ] && samtools faidx $ref
  
  args=($(xlsx2csv -i -d "\t" "${meta}" | tail -n +2 | awk '{print $1,$4,$10}' | egrep "${strains[$baseline]}" | awk '{print $3}' | sed "s/${baseline}//g"))
  samples=${args[@]/%/\|};samples=${samples// };samples=${samples::-1}

  R1_a=($(ls $reads | grep "R1" | grep -E $samples))
  R2_a=($(ls $reads | grep "R2" | grep -E $samples))
  
  N=5             #Number of analysis run in parallel
  (
  for ((i=0;i<${#R1_a[@]};i++));do    #loop through number of filenames in R1_a
    ((j=j%N)); ((j++==0)) && wait
    
##~~~~ singleSampleProcessing.sh: map paired-end reads to baseline reference and call variations ~~~~##    
    bash singleSampleProcessing.sh "${R1_a[i]}" "${R2_a[i]}" "${ref}" "${reads}" "${outdir}"&
  done
  wait
  )
  
done
wait

## check for truncation, sometimes with BAM files truncation can occur.
## BAM files where this happened and following up files have to be deleted and recreated from the start.
clear

proceed="yes"
for file in ${outdir}aligned_bwa/*/*.bam;do
  if [ ! $(samtools quickcheck -v $file | wc -l) == 0 ];then
    echo $file
    echo "remove ${file%/*}"
    #rm -r ${file%/*}
    proceed="no"
  fi
done

if [ ${proceed} == "no" ];then 
  echo "Remove the sample folders displayed above and rerun the pipeline."
  echo "Some BAM files are truncated, this can happen with big BAM files."
  echo "Remove manually or with 'rm -r [path to folder]'"
  echo "Exiting now."
  exit
fi
#===========================================================================================================#



#===========================================================================================================#
#                           Merge gVCFs into one multi-sample gVCF

for baseline in "${!strains[@]}";do
  echo ${strains[$baseline]}
  mkdir -p ${outdir}GATK_pipe_output/${strains[$baseline]}            #Create output directory
  
  if [ ! -s ${outdir}GATK_pipe_output/${strains[$baseline]}/cohort.g.vcf.gz ];then
    # get strain sample pool
    args=($(xlsx2csv -i --skipemptycolumns -d "\t" "${meta}" | tail -n +2 | awk '{print $1,$4,$10}' | egrep "${strains[$baseline]}" | awk '{print $3}'))
    sample_L=${args[@]/%/\|};sample_L=${sample_L// };sample_L=${sample_L::-1}
    
    # find all gVCF file paths of the strain, create input variable for CombineGCFs (-V sample1 -V sample2 -V sample3..N)
    merge_input=$(find . -name *.g.vcf.gz | egrep "${outdir}aligned_bwa" | egrep "${sample_L[*]}" | sort -t _ -k 3 | sed 's/\.\//-V /g')
     
##~~~~ GATK: CombineGVCFs ~~~~##    
    java -Xmx25g -XX:ParallelGCThreads=8 -jar tools/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar CombineGVCFs \
      -R "${outdir}/shovill/${strains[$baseline]}/contigs.fa" \
      -O ${outdir}GATK_pipe_output/${strains[$baseline]}/cohort.g.vcf \
      $merge_input

    bgzip ${outdir}GATK_pipe_output/${strains[$baseline]}/cohort.g.vcf                    #bgzip ../cohort.g.vcf
    tabix -f -p vcf ${outdir}GATK_pipe_output/${strains[$baseline]}/cohort.g.vcf.gz       #index cohort.g.vcf.gz
  fi
done
#===========================================================================================================#



#===========================================================================================================#
#                           Filter and annotate multiVCF

for baseline in "${!strains[@]}";do
  (
  strain=${strains[$baseline]}                                  #Strain name

  mkdir -p ${outdir}GATK_pipe_output/$strain                    #Create output directory
  ref="${outdir}/shovill/${strains[$baseline]}/contigs.fa"      #Path to consensus reference
  
  bash multiSampleProcessing.sh $strain $ref ${outdir}/GATK_pipe_output/$strain
  )&
done
wait

exit
#                           End of script
#===========================================================================================================#





# make contig file linear, so grep -A1 "contigXXXXX" gets you the sequence.
while read line;do
  if [ "${line:0:1}" == ">" ];then 
    echo -e "\n"$line 
  else echo $line | tr -d '\n';fi
done < genome_data/toolOut/shovill/KPN9884/contigs.fa > genome_data/toolOut/shovill/KPN9884/linear_contigs.fa


<<COMMENT 
# create de novo assembly to use in BLAST search for a template genome
T0=[CLINICAL ISOLATE BASELINE T=0, C=0]
T0="DTK_18"
R1=$(ls ${reads}*R1* | egrep $T0)
R2=$(echo ${R1} | sed 's/R1.filt.fastq.gz$/R2.filt.fastq.gz/g')
if [ ! -f assemblies/${T0}_deNovo/contigs.fasta ];then
 tools/SPAdes-3.13.0-Linux/bin/spades.py -1 ${R1} -2 ${R2} -o ${T0}_deNovo -t 32 -k 21,33,55,77 --careful
 #sed -e 's/\(^>.*$\)/#\1#/' ${assembly}${T0}_deNovo/contigs.fasta | tr -d "\r" | tr -d "\n" | sed -e's/$/#/' |\
 #tr "#" "\n" | sed -e '/^$/d' | egrep -A1 ">NODE_1_" > ${assembly}${T0}_deNovo/node1.fasta
fi

#NZ_CP014696.2 (DTK_35) https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP014696.2
#CP031795.1 (DTK_18) https://www.ncbi.nlm.nih.gov/nucleotide/CP031795.1
#CP034281.1 (DTK_01) https://www.ncbi.nlm.nih.gov/nucleotide/CP034281.1
COMMENT
#declare -A references=( ["DTK_01"]="CP034281.1" ["DTK_18"]="CP031795.1" )

  # download NCBI reference genome  #redudant, is in assembler
  if [ ! -f ${consensus}${baseline}_NCBItemplate.fa ];then
    esearch -db nucleotide -query ${references[$baseline]} | efetch -db nucleotide -format fasta > ${consensus}${baseline}_NCBItemplate.fa
  fi