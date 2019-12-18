#!/bash/bin

#===========================================================================================================#
#                           IMPORT VARIABLES & SET PATH
#===========================================================================================================#
# T0 SAMPLE, NCBI REFERENCE, FW_READS, RV_READS
T0=${1}
reference=${2}
R1=${3}
R2=${4}
strain=$5

# add to PATH for Prokka & snpEff
PATH=${PATH}:/home/linuxbrew/.linuxbrew/bin/

# Rename the numeric variables to their original names, for better script readabiliy.
Reads="genome_data/raw_sequences/"
consensus="genome_data/toolOut/consensus/"

# and make directories if they dont exist
mkdir -p ${consensus}



#===========================================================================================================#
#                           CREATE BASELINE REFERENCES
#===========================================================================================================#
# set NCBI reference variable
REF=${consensus}${T0}_NCBItemplate


##~~~~ Bowtie2/samtools/seqtk: index and build consensus.fa reference  ~~~~##
if [ ! -f ${consensus}${T0}_consensus.fa ];then
  # build bowtie2 database
  bowtie2-build $REF.fa $REF
  
  if [ ! -f ${consensus}${T0}_mapped.bam ];then
    # bowtie2 mapping
    bowtie2 -x $REF -1 $R1 -2 $R2 --very-sensitive -p 32 -S ${consensus}${T0}_mapped.sam
    # samtools: sort.sam and convert to .bam
    samtools view -bS ${consensus}${T0}_mapped.sam | samtools sort - -o ${consensus}${T0}_mapped.bam
  fi
  
  # create vcf file
  if [ ! -f ${consensus}${T0}.vcf ];then
    samtools mpileup -uf $REF.fa ${consensus}${T0}_mapped.bam > ${consensus}${T0}.vcf
  fi
  
  if [ ! -f ${consensus}${T0}_consensus.fastq ];then
    bcftools call -c --ploidy 1 ${consensus}${T0}.vcf | vcfutils.pl vcf2fq > ${consensus}${T0}_consensus.fastq
  fi
  
  # convert .fastq to .fasta and set bases of quality lower than 20 to N
  seqtk seq -aQ64 -q20 -n N ${consensus}${T0}_consensus.fastq > ${consensus}${T0}_consensus.fa
fi

##~~~~ prokka: annotate consensus.fa reference  ~~~~##
if [ ! -d genome_data/toolOut/prokka/$strain ];then
  prokka --force --addgenes --prefix $strain --outdir genome_data/toolOut/prokka/$strain ${consensus}${T0}_consensus.fa
fi


##~~~~ snpEff: build database for baseline consensus.fa references  ~~~~##
if [ ! -f /home/linuxbrew/.linuxbrew/share/snpeff/data/$strain/genes.gbk ];then
  # create needed directory
  mkdir -p /home/linuxbrew/.linuxbrew/share/snpeff/data/$strain
  
  # copy genbank file and gff file for creating the database
  cp genome_data/toolOut/prokka/$strain/$strain.gbk /home/linuxbrew/.linuxbrew/share/snpeff/data/$strain
  mv /home/linuxbrew/.linuxbrew/share/snpeff/data/$strain/$strain.gbk /home/linuxbrew/.linuxbrew/share/snpeff/data/$strain/genes.gbk
  
  cp genome_data/toolOut/prokka/$strain/$strain.gff /home/linuxbrew/.linuxbrew/share/snpeff/data/$strain
  mv /home/linuxbrew/.linuxbrew/share/snpeff/data/$strain/$strain.gff /home/linuxbrew/.linuxbrew/share/snpeff/data/$strain/regulation.gff
  
  # add database to config file, if it isn't in there yet
  if ! grep -Fxq "$strain.genome : $strain" /home/linuxbrew/.linuxbrew/share/snpeff/snpEff.config;then
    echo "$strain.genome : $strain" >> /home/linuxbrew/.linuxbrew/share/snpeff/snpEff.config
  fi
  
  # build snpEff database
  snpEff build -genbank -v $strain
fi

exit

