#!/bin/bash
#hostname
module load blast/blast-2.6.0
module load gcc/gcc-7.2.0 python/python-anaconda3.5 
module load racon
module load samtools/samtools-1.3.3
module load spades/spades-3.11.0
module load java
module load bowtie2/bowtie2-2.3.4 
module load bedtools/bedtools2.25

phix_db_path=/urigo/DBs/phix/PhiX
SEQDIR_PE=/urigo/4N_10N_16N_24N_48N/
SEQDIR_Pac=/scratch200/urineri/Not_my_raw_reads/pacbio/pacbio/
OUTDIR=/scratch200/urineri/hybrid/
sample_List=(Sample_48N Sample_4N)
bbduk="/scratch200/urineri/resources/bbmap/bbduk.sh"
THREADS=8 #use more
contam=/scratch200/urineri/ccSAG/example/known_cont.fasta

#Fastqc Raw Reads
for sample in ${sample_List[*]};
do 
  sample_name=$(echo $sample | cut -d'_' -f 2) #remove the "Sample" 
  echo "Generating intial FastQC report for: " $sample_name
  mkdir $OUTDIR/$sample_name
  curr_out=$OUTDIR/$sample_name
  mkdir $curr_out/FastQC_raw_reads/
  cd $curr_out/FastQC_raw_reads/
    fastqc -t $THREADS $SEQDIR_PE/$sample/"$sample_name"_R1.fastq -o $curr_out/FastQC_raw_reads/
    fastqc -t $THREADS $SEQDIR_PE/$sample/"$sample_name"_R2.fastq -o $curr_out/FastQC_raw_reads/
    echo "Finshed intial QC report for: " $sample_name
cd $SEQDIR_PE
done

for sample in ${sample_List[*]};
do 
  sample_name=$(echo $sample | cut -d'_' -f 2) #remove the "Sample" 
  echo "Started mapping to phiX174 for: " $sample_name
  #mkdir $OUTDIR/$sample_name
  curr_out=$OUTDIR/$sample_name
  cd $curr_out
    bowtie2 -x $phix_db_path -1 $SEQDIR_PE/$sample/"$sample_name"_R1.fastq  -2 "$SEQDIR_PE"/$sample/"$sample_name"_R2.fastq  -S $curr_out/"$sample_name"_to_phix174.sam --threads $THREADS
  samtools view -@ $THREADS -bS $curr_out/"$sample_name"_to_phix174.sam > $curr_out/"$sample_name"_to_phix174_mapped_and_unmapped.bam
  samtools view -b -@ $THREADS -f 12 -F 256 $curr_out/"$sample_name"_to_phix174_mapped_and_unmapped.bam > $curr_out/"$sample_name"_to_phix174_both_EndsUnmapped.bam
  samtools sort -@ $THREADS -n $curr_out/"$sample_name"_to_phix174_both_EndsUnmapped.bam -o $curr_out/"$sample_name"_to_phix174_both_EndsUnmapped_sorted.bam
  bedtools bamtofastq -i  $curr_out/"$sample_name"_to_phix174_both_EndsUnmapped_sorted.bam -fq $curr_out/"$sample_name"_R1.fastq -fq2 $curr_out/"$sample_name"_R2.fastq
  echo "Finshed: " $sample_name
  rm ./*.bam
  rm ./*.sam
cd $SEQDIR_PE
done

#Trimming PE
for sample in ${sample_List[*]};
do 
  sample_name=$(echo $sample | cut -d'_' -f 2) #remove the "Sample" 
  echo "Started trimmomatic for: " $sample_name
  #mkdir $OUTDIR/$sample_name
  curr_out=$OUTDIR/$sample_name
  cd $curr_out
  java -jar /powerapps/sources/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads $THREADS -phred33 -trimlog $curr_out/logfile $curr_out/"$sample_name"_R1.fastq $curr_out/"$sample_name"_R2.fastq $curr_out/"$sample_name"_trim_paired_R1.fastq $curr_out/"$sample_name"_trim_unpaired_R1.fastq $curr_out/"$sample_name"_trim_paired_R2.fastq $curr_out/"$sample_name"_trim_unpaired_R2.fastq ILLUMINACLIP:/powerapps/sources/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  echo "Finshed trimmomatic for: " $sample_name
cd $SEQDIR_PE
done

#Trimming PacBio (make sure to merge the input subreads.fastq to one file)
for sample in ${sample_List[*]};
do 
  sample_name=$(echo $sample | cut -d'_' -f 2) #remove the "Sample" 
  echo "started trimming PacBio for: " $sample_name
  curr_out=$OUTDIR/$sample_name
  cd $curr_out
  $bbduk in=$SEQDIR_Pac/"$sample_name"/"$sample_name"_cated_pacbio.fastq out=$curr_out/"$sample_name"_filtered.fq minlen=100 maq=8 threads=$THREADS
  echo "Finshed trimming PacBio: " $sample_name
  cd $SEQDIR_Pac
done

#Fastqc post trimmomatic Illumina reads (as well as post bbduk PacBio reads). 
for sample in ${sample_List[*]};
do 
  sample_name=$(echo $sample | cut -d'_' -f 2) #remove the "Sample" 
  echo "Generating FastQC report for reads post trimmomatic/bbduk: " $sample_name
  #mkdir $OUTDIR/$sample_name
  curr_out=$OUTDIR/$sample_name
  mkdir $curr_out/FastQC_post_trim_reads/
  cd $curr_out/FastQC_post_trim_reads/
    fastqc -t $THREADS $curr_out/"$sample_name"_trim_paired_R1.fastq -o $curr_out/FastQC_post_trim_reads/
    fastqc -t $THREADS $curr_out/"$sample_name"_trim_paired_R2.fastq -o $curr_out/FastQC_post_trim_reads/
    fastqc -t $THREADS $curr_out/"$sample_name"_filtered.fq -o $curr_out/FastQC_post_trim_reads/
    echo "Finshed post trimming QC report for: " $sample_name
cd $SEQDIR_PE
done

#Assembler (unicycler, consider using Bold mode instead of default (Normal mode).
for sample in ${sample_List[*]};
do 
  sample_name=$(echo $sample | cut -d'_' -f 2) #remove the "Sample" 
  echo "Started assembly for: " $sample_name
  curr_out=$OUTDIR/$sample_name
  cd $curr_out
  mkdir ./unicyclr
  unicycler -1 $curr_out/"$sample_name"_trim_paired_R1.fastq -2 $curr_out/"$sample_name"_trim_paired_R2.fastq -l $curr_out/*_filtered.fq --threads $THREADS --keep 3 --verbosity 2 --pilon_path /scratch200/urineri/pilon-1.23.jar -o ./unicyclr/

  echo "Finshed: " $sample_name
cd $OUTDIR
done
