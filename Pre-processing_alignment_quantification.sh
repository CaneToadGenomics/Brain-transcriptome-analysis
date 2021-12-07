## Check reads quality (FastQC 0.11.7)
for i in $(cat /path/samples.txt)
do
fastqc /path/${i}_1.fastq.gz --threads 2 --extract --outdir=/path ${i}_1.fastq
fastqc /path/${i}_2.fastq.gz --threads 2 --extract --outdir=/path ${i}_2.fastq
done

## Trim adapters and low-quality reads (Trimmomatic 0.38). Change parameters based on FastQC outcome.
for i in $(cat /path/samples.txt)
do
java -jar /apps/trimmomatic/0.38/trimmomatic-0.38.jar PE /path/${i}_1.fastq.gz /path/${i}_2.fastq.gz /path/${i}_1_paired.fastq.gz /path/${i}_1_unpaired.fastq.gz /path/${i}_2_paired.fastq.gz /path/${i}_2_unpaired.fastq.gz ILLUMINACLIP:/apps/trimmomatic/0.38/adapters/TruSeq3-PE.fa:2:30:10:4 HEADCROP:13 AVGQUAL:30 MINLEN:36
done

## Check reads quality again (FastQC 0.11.7)
for i in $(cat /path/samples.txt)
do
fastqc /path/${i}_1_paired.fastq.gz --threads 2 --extract --outdir=/path ${i}_1_paired.fastq
fastqc /path/${i}_2_paired.fastq.gz --threads 2 --extract --outdir=/path ${i}_2_paired.fastq
done

## Build genome index (Star 2.7.2b). Change parameters accordingly.
STAR --runThreadN 2 —-runMode genomeGenerate --genomeDir /path/genome -- genomeFastaFiles /path/Transcriptome.fa --sjdbGTFfile /path/Transcriptome.gff --sjdbOverhang 112

## Map reads to transcriptome (Star 2.7.2b).
for i in $(cat /path/samples.txt)
do
STAR --runThreadN 2 --genomeDir /path/genome --readFilesIn /path/${i}_1_paired.fastq /path/${i}_2_paired.fastq --outFileNamePrefix /path/${i} --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 2 --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic
done

## Quantify gene expression (Salmon 1.2.1).
for i in $(cat /path/samples.txt)
do
salmon quant -t /path/Transcriptome.fa -l A -a /path/${i}Aligned.sortedByCoord.out.bam -o /path/${i}_quant
done



