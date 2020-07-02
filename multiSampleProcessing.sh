#!/bash/bin
#
#===========================================================================================================#
#                           Introduction
usage(){
clear
echo "
multiSampleProcessing version 0.0.1 part of the runNGSdownstreamPipeline (build 2020-02-18), by Yob Haakman <yob1997@live.nl>.

Usage:  bash multiSampleProcessing.sh [STRAIN] [PATH_TO_REFERENCE_GENOME.fa] [OUTPUT_DIRECTORY]

Required parameters:
  [STRAIN]                      :Strain name, will be used in annotation and for building a snpEff database
  [PATH_TO_REFERENCE_GENOME.fa] :Path to reference genome (fasta format)
  [OUTPUT_DIRECTORY]            :Folder to output results to [stdout]

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

# add to PATH for VT and snpEff
PATH=${PATH}:/home/linuxbrew/.linuxbrew/bin/

strain=$1
ref=$2
outdir=$3
#===========================================================================================================#



#===========================================================================================================#
#                           SNP calling with GATK

if [ ! -s ${outdir}/raw.vcf.gz ];then
  java -Xmx25g -XX:ParallelGCThreads=32 -jar tools/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar GenotypeGVCFs \
    -R ${ref} \
    -ploidy 1 \
    --standard-min-confidence-threshold-for-calling 0 \
    -V ${outdir}/cohort.g.vcf.gz \
    -O ${outdir}/raw.vcf.gz
fi
#===========================================================================================================#



#===========================================================================================================#
#                           Decompose & normalize SNPs with VT

if [ ! -s ${outdir}/decompose.vcf.gz ];then
  vt decompose -s -o ${outdir}/decompose.vcf.gz ${outdir}/raw.vcf.gz
  tabix -f -p vcf ${outdir}/decompose.vcf.gz
fi

if [ ! -s ${outdir}/normalize.vcf.gz ];then
  vt normalize -r ${ref} \
    -o ${outdir}/normalize.vcf.gz ${outdir}/decompose.vcf.gz
  tabix -f -p vcf ${outdir}/normalize.vcf.gz
fi
#===========================================================================================================#



#===========================================================================================================#
#                           SNP filtering with GATK

##~~~~ GATK: SelectVariants extract SNPs ~~~~##
if [ ! -s ${outdir}/normalized_snps.vcf.gz ];then
  java -Xmx25g -XX:ParallelGCThreads=32 -jar tools/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar SelectVariants \
    -R ${ref} \
    -V ${outdir}/normalize.vcf.gz \
    -O ${outdir}/normalized_snps.vcf.gz \
    -select-type SNP
fi

##~~~~ GATK: VariantFiltrations filter SNPs ~~~~##
# (https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set)
##QD, Description="Variant Confidence/Quality by Depth">
##FS, Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##MQ, Description="RMS Mapping Quality">
##MQRankSum, Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##ReadPosRankSum, Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##SOR, Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
if [ ! -s ${outdir}/filter_snps.vcf.gz ];then
  java -Xmx25g -XX:ParallelGCThreads=32 -jar tools/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar VariantFiltration \
    -R ${ref} \
    -V ${outdir}/normalized_snps.vcf.gz \
    -O ${outdir}/filter_snps.vcf.gz \
    --filter-name "SNP-QD2" -filter-expression "QD < 2.0" \
    --filter-name "SNP-FS60" -filter-expression "FS > 60.0" \
    --filter-name "SNP-MQ40" -filter-expression "MQ < 40.0" \
    --filter-name "SNP-MQRankSum-12.5" -filter-expression "MQRankSum > -12.5" \
    --filter-name "SNP-ReadPosRankSum-8" -filter-expression "ReadPosRankSum > -8.0" \
    --filter-name "SNP-SOR3" -filter-expression "SOR > 3.0"
    
fi

##~~~~ GATK: SelectVariants extract INDELs ~~~~##
if [ ! -s ${outdir}/normalized_indels.vcf.gz ];then
  java -Xmx25g -XX:ParallelGCThreads=32 -jar tools/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar SelectVariants \
    -R ${ref} \
    -V ${outdir}/normalize.vcf.gz \
    -O ${outdir}/normalized_indels.vcf.gz \
    -select-type INDEL
fi

##~~~~ GATK: VariantFiltrations filter INDELs ~~~~##
# (https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set)
##QD, Description="Variant Confidence/Quality by Depth">
##FS, Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##ReadPosRankSum, Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##SOR, Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
if [ ! -s ${outdir}/filter_indels.vcf.gz ];then
  java -Xmx25g -XX:ParallelGCThreads=32 -jar tools/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar VariantFiltration \
    -R ${ref} \
    -V ${outdir}/normalized_indels.vcf.gz \
    -O ${outdir}/filter_indels.vcf.gz \
    --filter-name "INDEL-QD2" -filter-expression "QD < 2.0" \
    --filter-name "INDEL-FS200" -filter-expression "FS > 200.0" \
    --filter-name "INDEL-ReadPosRankSum-20" -filter-expression "ReadPosRankSum > -20.0" \
    --filter-name "INDEL-SOR10" -filter-expression "SOR > 10.0"

fi

##~~~~ GATK: MergeVcfs put all variants back into one VCF ~~~~##
if [ ! -s ${outdir}/merged_filt.vcf.gz ];then
  java -Xmx25g -XX:ParallelGCThreads=32 -jar tools/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar MergeVcfs \
    -I ${outdir}/filter_snps.vcf.gz \
    -I ${outdir}/filter_indels.vcf.gz \
    -O ${outdir}/merged_filt.vcf.gz
fi
#===========================================================================================================#




#~~~~ snpEff: annotate VCF with snpEff and Snpsift ~~~~##
if [ ! -s  ${outdir}/annotated.vcf.gz ];then
  set -o pipefail
  # annotate and manipulate with snpEff and SnpSift
  snpEff -no-downstream -no-upstream \
  ${strain}_pgap ${outdir}/merged_filt.vcf.gz > /dev/stdout |\
  SnpSift varType /dev/stdin > ${outdir}/annotated.vcf

  # move snpEff genes.txt and summary.html 
  mv snpEff_genes.txt ${outdir}/snpEff_genes.txt
  mv snpEff_summary.html ${outdir}/snpEff_summary.html

  # bgzip annotated file and create index
  bgzip ${outdir}/annotated.vcf
  tabix -f -p vcf ${outdir}/annotated.vcf.gz
fi
exit
#                           End of script
#===========================================================================================================#