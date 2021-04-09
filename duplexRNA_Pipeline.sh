#!/bin/bash
echo "############################"
echo -e "Duplex RNASeq pipeline version v1.0 initiated.\n"

if [ $# -ne 2 ]
then
        echo "Usage: <Project Folder> <sampleID>"
       exit 65
fi

# Data directory and sample
data_dir=$1
sample=$2
pipeline_dir=$3
resources_dir=$4

echo $data_dir
echo $sample

#Directories
fullpath=`readlink -e "${BASH_SOURCE[0]}"`
echo -e "Full path is $fullpath"
dirn=`dirname $fullpath`
parent_dir=`dirname $dirn`
echo -e "fullpath=$fullpath\tdirn=$dirn\tparent dir=$parent_dir"

fastq_dir=$data_dir/fastq
qc_dir=$data_dir/qc
fastqc_dir=$data_dir/qc/fastqc

mkdir -p $qc_dir
mkdir -p $fastqc_dir

#Reference files
genome_file=${resources_dir}/PlasmoDB-34_Pfalciparum3D7_Genome.fasta
star_genome_folder=${resources_dir}/Plasmodium_v34_STAR

echo "############################"
time_start_english=`date`
master_time_start=`date +%s`
time_start=`date +%s`
echo -e "SAMPLE\t$sample"
echo "$sample duplex RNASeq pipeline v1.0 initiated at $time_start_english"
echo "############################"
time_end=`date +%s`
(time_exec=`expr $(( $time_end - $time_start ))`; echo "$sample read count  completed in $time_exec seconds\n")
time_start=`date +%s`

echo "############################"
echo -e "Running fastp\n"

# fastp path
fastp -i $fastq_dir/${sample}_R1_001.fastq.gz -I $fastq_dir/${sample}_R2_001.fastq.gz -o $fastq_dir/${sample}_R1_trimmed.fastq.gz -O $fastq_dir/${sample}_R2_trimmed.fastq.gz -h $qc_dir/${sample}_fastp.html -l 15

echo "############################"
time_end=`date +%s`
(time_exec=`expr $(( $time_end - $time_start ))`; echo "$sample Trimming completed in $time_exec seconds")
time_start=`date +%s`

echo "############################"
echo -e "Run Fastqc on the raw  reads"

fastqc $fastq_dir/$sample"_R1_trimmed.fastq.gz" -o $fastqc_dir --extract                                                                                                  
fastqc $fastq_dir/$sample"_R2_trimmed.fastq.gz" -o $fastqc_dir --extract                                                                                           

echo "############################"
time_end=`date +%s`
(time_exec=`expr $(( $time_end - $time_start ))`; echo "$sample FastQC on fastq completed in $time_exec seconds")
time_start=`date +%s`

echo "############################"
echo -e "Run STAR"

star_dir=$data_dir/STAR
mkdir -p $star_dir

STAR --genomeDir $star_genome_folder \
    --runThreadN 16 \
    --readFilesIn $fastq_dir/$sample"_R1_trimmed.fastq.gz" $fastq_dir/$sample"_R2_trimmed.fastq.gz" \
    --outFileNamePrefix $star_dir/$sample.star. \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbGTFfile ${resources_dir}/PlasmoDB-34_Pfalciparum3D7_ForSTAR.gtf \
   --genomeSAindexNbases 11 \
   --alignIntronMax 500

echo "############################"
time_end=`date +%s`
(time_exec=`expr $(( $time_end - $time_start ))`; echo "$sample STAR Alignments completed in $time_exec seconds")
time_start=`date +%s`

echo "############################"
echo -e "Index the file"

star_sorted_file=$star_dir/$sample.star.Aligned.sortedByCoord.out.bam

samtools index $star_sorted_file
java -Xmx1g -jar $PICARD SortSam SO=coordinate INPUT=$star_dir/$sample".starAligned.out.sam" OUTPUT=$bam_dir/$sample".star.sorted.bam" CREATE_INDEX=true

echo "############################"
time_end=`date +%s`
(time_exec=`expr $(( $time_end - $time_start ))`; echo "$sample STAR BAM Indexing  completed in $time_exec seconds")
time_start=`date +%s`

echo "############################"
echo -e "IGVTools toTDF"

tdf_dir=$data_dir/TDF
mkdir -p $tdf_dir
${IGVTOOLS}/igvtools count -z 5 -w 5 --strands first $star_sorted_file $tdf_dir/$sample".star.sorted.stranded.tdf" ${genome_file}
${IGVTOOLS}/igvtools count -z 5 -w 5 --strands first --minMapQuality 5 $star_sorted_file $tdf_dir/$sample".star.sorted.stranded.mapq5.tdf" ${genome_file}

echo "############################"
time_end=`date +%s`
(time_exec=`expr $(( $time_end - $time_start ))`; echo "$sample toTDF completed in $time_exec seconds")
time_start=`date +%s`

echo "###############################"
echo -e "Filter low quality reads out"

python $pipeline_dir/remove_bad_reads.py --input $star_sorted_file --out ${star_sorted_file}_filtered.bam
samtools index ${star_sorted_file}_filtered.bam

echo "############################"
time_end=`date +%s`
(time_exec=`expr $(( $time_end - $time_start ))`; echo "$sample Filtering the BAM completed in $time_exec seconds")
time_start=`date +%s`

echo "############################"
echo -e "IGVTools toTDF"

tdf_dir=$data_dir/TDF
mkdir -p $tdf_dir
IGVTOOLS=/projects/b1103/software/IGVTools_2.4.19 
${IGVTOOLS}/igvtools count -z 5 -w 5 --strands first ${star_sorted_file}_filtered.bam $tdf_dir/$sample".star.filtered.sorted.stranded.tdf" ${genome_file}          

echo "############################"
time_end=`date +%s`
(time_exec=`expr $(( $time_end - $time_start ))`; echo "$sample Filtered BAM toTDF completed in $time_exec seconds")
time_start=`date +%s`

echo "############################"
echo -e "Run Picard Insert size distribution"

insert_dir=${qc_dir}/insert_sizes
mkdir -p $insert_dir           
java -Xmx5g -jar $PICARD CollectInsertSizeMetrics HISTOGRAM_FILE=$insert_dir/$sample".star.sorted.insertsizes.histogram.pdf" INPUT=$star_sorted_file OUTPUT=$insert_dir/$sample"star.sorted.insertsizes.txt" 
java -Xmx5g -jar $PICARD CollectInsertSizeMetrics HISTOGRAM_FILE=$insert_dir/$sample".star.filtered.sorted.insertsizes.histogram.pdf" INPUT=${star_sorted_file}_filtered.bam OUTPUT=$insert_dir/$sample"star.filtered.sorted.insertsizes.txt"

echo "############################"
time_end=`date +%s`
(time_exec=`expr $(( $time_end - $time_start ))`; echo "$sample Insert size metrics completed in $time_exec seconds")
time_start=`date +%s`

echo "###########################"
echo -e " RUNNING FEATURECOUNTS\n"

featurecount_dir=$data_dir/FeatureCount
mkdir -p $featurecount_dir

# JULY 10TH 2020
# ADDING NEW GTF FILE WITH NEWER BETTER ANNOTATIONS
numThreads=16

# -p for paired end
# -O to assign reads to overlapping features
# -t gene to output reads at gene level
# -s 1 because strandedness is in first of pair read (checked in IGV)
# -C does not count reads in which one mate is mapped to one chromosome and other to another chromosome
featureCounts -p -O --largestOverlap -C -t gene -s 1 -a $refgtf -g ID -o $featurecount_dir/${sample}.filtered.featurecounts.txt -T $numThreads ${star_sorted_file}_filtered.bam

echo "############################"
time_end=`date +%s`
(time_exec=`expr $(( $time_end - $time_start ))`; echo "$sample FeatureCounts completed in $time_exec seconds")
time_start=`date +%s`

