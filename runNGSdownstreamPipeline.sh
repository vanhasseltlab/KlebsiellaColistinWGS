#!/bash/bin
#ready for git!

#read -p "Set output directory" outdir
outdir="genome_data/toolOut/"

#===========================================================================================================#
#                           Create baseline references
#===========================================================================================================#
reads="genome_data/raw_sequences/"            # Path to filtered sample reads
consensus="${outdir}consensus/"    # Path to consensus sequences used as reference


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$    PROJECT SPECIFIC    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
declare -A references=( ["DTK_01"]="CP034281.1" ["DTK_18"]="CP031795.1" ["DTK_35"]="NZ_CP014696.2" )
declare -A strains=(["DTK_01"]="KPN9884" ["DTK_18"]="KPN9749" ["DTK_35"]="KPN700603")

# create de novo assembly to use in BLAST search for a template genome
#T0=[CLINICAL ISOLATE BASELINE T=0, C=0]
#R1=$(ls ${Reads}*R1* | egrep $T0)
#R2=$(echo ${R1} | sed 's/R1.filt.fastq.gz$/R2.filt.fastq.gz/g')
#if [ ! -f ${assembly}${T0}_deNovo/contigs.fasta ];then
# tools/SPAdes-3.13.0-Linux/bin/spades.py -1 ${R1} -2 ${R2} -o ${T0}_deNovo -t 32 -k 21,33,55,77 --careful
# sed -e 's/\(^>.*$\)/#\1#/' ${assembly}${T0}_deNovo/contigs.fasta | tr -d "\r" | tr -d "\n" | sed -e's/$/#/' |\
# tr "#" "\n" | sed -e '/^$/d' | egrep -A1 ">NODE_1_" > ${assembly}${T0}_deNovo/node1.fasta
#fi

#NZ_CP014696.2 (DTK_35) https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP014696.2
#CP031795.1 (DTK_18) https://www.ncbi.nlm.nih.gov/nucleotide/CP031795.1
#CP034281.1 (DTK_01) https://www.ncbi.nlm.nih.gov/nucleotide/CP034281.1

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$    PROJECT SPECIFIC    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


for baseline in "${!references[@]}";do
  # download NCBI reference genome
  if [ ! -f ${consensus}${baseline}_NCBItemplate.fa ];then
    esearch -db nucleotide -query ${references[$baseline]} | efetch -db nucleotide -format fasta > ${consensus}${baseline}_NCBItemplate.fa
  fi
  
  # create consensus sequences of baseline samples
  if [ ! -f ${consensus}${baseline}_consensus.fa ];then
    R1=$(ls ${reads}*R1* | egrep ${baseline})
    R2=$(echo ${R1} | sed 's/R1.filt.fastq.gz$/R2.filt.fastq.gz/g')
    bash createBaseLineREF.sh $baseline ${references[$baseline]} $R1 $R2 ${strains[$baseline]}
  fi
done


clear
echo "Baseline references created!"
echo "Next steps can take a long (-+ 2 days for 200 samples) time to run."

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
#                           Mapping to baseline references
#===========================================================================================================#
for baseline in "${!references[@]}";do
  echo $baseline

  R1=$(ls ${reads}*R1* | egrep ${baseline})                             #Path to forward reads
  R2=$(echo ${R1} | sed 's/R1.filt.fastq.gz$/R2.filt.fastq.gz/g')       #Path to reverse reads
  ref="${outdir}consensus/${baseline}_consensus.fa"                     #Path to reference consensus genome
  
##~~~~ BWA MEM/GATK/samtools: index reference genome ~~~~##
  if [ ! -s $ref.ann ];then bwa index $ref;fi
  if [ ! -s ${outdir}consensus/${baseline}_consensus.dict ];then java -jar tools/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar CreateSequenceDictionary -R $ref;fi
  if [ ! -s $ref.fai ];then samtools faidx $ref;fi
  

  args=($(xlsx2csv -i -d "\t" 'Sample overview_seq.xlsx' | tail -n +2 | awk '{print $1,$4,$10}' | egrep "${strains[$baseline]}" | awk '{print $3}'))
  samples=${args[@]/%/\|};samples=${samples// };samples=${samples::-1}


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$    PROJECT SPECIFIC    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  # add clinical samples to strain sample pool
  if [ ${strains[$baseline]} == "KPN9749" ];then samples=$(echo $samples"|FSTK_120|FSTK_121");fi
  if [ ${strains[$baseline]} == "KPN9884" ];then samples=$(echo $samples"|FSTK_122|FSTK_123");fi
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$    PROJECT SPECIFIC    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


  R1_a=($(ls $reads | grep "R1" | grep -E $samples | tail -n +2))
  R2_a=($(ls $reads | grep "R2" | grep -E $samples | tail -n +2))
  
  N=5             #Number of analysis run in parallel
  (
  for ((i=0;i<${#R1_a[@]};i++));do    #loop through number of filenames in R1_a
    ((j=j%N)); ((j++==0)) && wait
    
##~~~~ singleSampleProcessing.sh: map paired-end reads to baseline reference and call variations ~~~~##    
    bash singleSampleProcessing.sh "${R1_a[i]}" "${R2_a[i]}" "${ref}" "${reads}"&
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
  echo "Exiting now."
  exit
fi



#===========================================================================================================#
#                           Merge gVCFs into one multi-sample gVCF
#===========================================================================================================#
for baseline in "${!strains[@]}";do
  echo ${strains[$baseline]}
  mkdir -p ${outdir}GATK_pipe_output/${strains[$baseline]}            #Create output directory
  
  if [ ! -s ${outdir}GATK_pipe_output/${strains[$baseline]}/cohort.g.vcf.gz ];then
    # get strain sample pool
    args=($(xlsx2csv -i -d "\t" 'Sample overview_seq.xlsx' | tail -n +2 | awk '{print $1,$4,$10}' | egrep "${strains[$baseline]}" | awk '{print $3}'))
    sample_L=${args[@]/%/\|};sample_L=${sample_L// };sample_L=${sample_L::-1}
    
    # add clinical samples to strain sample pool
    if [ ${strains[$baseline]} == "KPN9749" ];then sample_L=$(echo $sample_L"|FSTK_120|FSTK_121");fi
    if [ ${strains[$baseline]} == "KPN9884" ];then sample_L=$(echo $sample_L"|FSTK_122|FSTK_123");fi
    
    # find all gVCF file paths of the strain, create input variable for CombineGCFs (-V sample1 -V sample2 -V sample3..N)
    merge_input=$(find . -name *.g.vcf.gz | grep "toolOut" | grep -E "${sample_L[*]}" | sort -t _ -k 3 | sed 's/\.\//-V /g')
    
##~~~~ GATK: CombineGVCFs ~~~~##    
    java -Xmx25g -XX:ParallelGCThreads=8 -jar tools/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar CombineGVCFs \
      -R ${outdir}consensus/${baseline}_consensus.fa \
      -O ${outdir}GATK_pipe_output/${strains[$baseline]}/cohort.g.vcf \
      $merge_input

    bgzip ${outdir}GATK_pipe_output/${strains[$baseline]}/cohort.g.vcf                    #bgzip ../cohort.g.vcf
    tabix -f -p vcf ${outdir}GATK_pipe_output/${strains[$baseline]}/cohort.g.vcf.gz       #index cohort.g.vcf.gz
  fi
  
  
done



#===========================================================================================================#
#                                     filter and annotate multiVCF
#===========================================================================================================#
for baseline in "${!strains[@]}";do
  (
  strain=${strains[$baseline]}                          #Strain name

  mkdir -p ${outdir}GATK_pipe_output/$strain            #Create output directory
  ref="${outdir}consensus/${baseline}_consensus.fa"     #Path to consensus reference
  
  bash multiSampleProcessing.sh $strain $ref ${outdir}GATK_pipe_output/$strain
  )&
done
wait



exit


echo "give path to reference genome"
reference=$(read -p)

