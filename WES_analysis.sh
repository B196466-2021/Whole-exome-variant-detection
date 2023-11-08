mkdir {raw,clean,align,genome,hg19_VCF}
# Run in the main project directory and generate two filtered fq files in the clean directory.
#Raw data quality control, using fastp
nohup fastp -i raw/wes.1.fq.gz -o clean/wes.1.clean.fq.gz -I raw/wes.2.fq.gz -O clean/wes.2.clean.fq.gz &
##Download the reference genome (hg19) file, store it in the genome directory, and create an index:
for i in $(seq 1 22) X Y M;
do 
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${i}.fa.gz
done

for i in $(seq 1 22) X Y M;
do cat chr${i}.fa.gz >> hg19.fa;
done

bwa index -a bwtsw -p hg19 hg19.fa

#alignment
# Run this command in the clean directory to generate the sam file in the align directory.
bwa mem -t 4 -M -R "@RG\tID:lane1\tPL:illumina\tLB:library\tSM:wes" /genome/hg19.fa  /raw/wes.1.fq  /raw/wes.2.fq > wes.sam
#Easy to store
samtools view -b -S wes.sam > wes.bam
#Sequential sorting to facilitate subsequent operations
samtools sort wes.bam -o wes.sorted.bam
#Statistical comparison information
samtools flagstat wes.sorted.bam > wes.sorted.bam.flagstat


##Exon region coverage
gatk CreateSequenceDictionary -R hg19.fa -O hg19.dict
gatk BedToIntervalList -I S31285117_Regions.bed -O Exon.Interval.bed -SD ./genome/hg19.dict

##Mark PCR repeats and index them
gatk --java-options "-Xmx10G -Djava.io.tmpdir=./" MarkDuplicates -I wes.sorted.bam -O wes.sorted.MarkDuplicates.bam -M wes.sorted.bam.metrics > log.mark 2>&1
samtools view -f 1024 wes.sorted.MarkDuplicates.bam | less
samtools index wes.sorted.MarkDuplicates.bam

#Variation detection
#Recalibrate base quality values (BQSR)
#First download the variant annotation file
wget  -c -r -nd -np -k -L -p ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19
#The downloaded files are all compressed and need to be decompressed before use
samtools faidx hg19.fa
GENOME=./genome/hg19.fa
hg19_VCF=./hg19_VCF/
gatk --java-options "-Xmx10G -Djava.io.tmpdir=./" BaseRecalibrator \
-R $GENOME -I ./align/wes.sorted.MarkDuplicates.bam \
--known-sites $hg19_VCF/1000G_phase1.indels.hg19.sites.vcf \
--known-sites $hg19_VCF/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
--known-sites $hg19_VCF/dbsnp_138.hg19.vcf \
-L S31285117_Regions.bed -O wes.recal_data.table


gatk --java-options "-Xmx10G -Djava.io.tmpdir=./" ApplyBQSR \
-R $GENOME -I ./align/wes.sorted.MarkDuplicates.bam \
-bqsr wes.recal_data.table -L S31285117_Regions.bed  \
-O wes.sorted.MarkDuplicates.BQSR.bam


# gatk AnalyzeCovariates -bqsr wes.recal_data.table -plots wes.recal_data.table.plot

##Mutation detection
# 1.Generate intermediate file gvcf
gatk --java-options "-Xmx8G -Djava.io.tmpdir=./" HaplotypeCaller -R $GENOME \
--emit-ref-confidence GVCF -I wes.sorted.MarkDuplicates.BQSR.bam \
-D $hg19_VCF/dbsnp_138.hg19.vcf -L S31285117_Regions.bed  -O wes.gvcf


# 2.Detecting mutations through GVCF
gatk --java-options "-Xmx8G -Djava.io.tmpdir=./" GenotypeGVCFs \
-R $GENOME -V wes.gvcf -L S31285117_Regions.bed  \
-O wes.raw.vcf

#Variation Quality Control and Filtering
gatk --java-options "-Xmx8G -Djava.io.tmpdir=./" VariantRecalibrator -R $GENOME -V wes.raw.vcf \
-resource hapmap,known=false,training=true,truth=true,prior=15.0:$hg19_VCF/hapmap_3.3.hg19.sites.vcf \
-resource omini,known=false,training=true,truth=false,prior=12.0:$hg19_VCF/1000G_omni2.5.hg19.sites.vcf \ 
-resource 1000G,known=false,training=true,truth=false,prior=10.0:$hg19_VCF/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
-resource dbsnp,known=true,training=false,truth=false,prior=6.0:$hg19_VCF/dbsnp_138.hg19.vcf \
-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode SNP \
-O wes.snps.recal.vcf --tranches-file wes.snps.tranches --rscript-file wes.snps.plots.R & 


gatk  ApplyRecalibration  -V wes.raw.vcf -O wes.VQSR.vcf --recal-file wes.snps.recal.vcf  \
--tranches-file wes.snps.tranches  -mode SNP

##Mutation annotation
#Database download
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
 # -buildver 
 # -downdb 
 # -webfrom annovar 
 # humandb/ 

 #SNP annotation
 table_annovar.pl snp.avinput $humandb -buildver hg19 -out snpanno \
-remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all \
-operation g,r,r,f -nastring . -csvout