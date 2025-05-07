# login
ssh group_2@vpsmbw-p-srv04.mbw.su.se # password: group_2
# back to home dir
cd ~
# our sequencing data
/home/admin4/2025_yeast_NGS/group_2_S2_L001_R1_001.fastq.gz


# copy ref genome, sample and control into our folder
cd /home/group_2/uz_real_data/1_raw_data
cp /home/admin4/2025_yeast_NGS/group_2_S2_L001_R1_001.fastq.gz .
cp /home/admin4/2025_yeast_NGS/control1_R1.fastq .
cp /home/admin4/SNP_calling/first_part/reference/S288C_reference_plus_2u.fa .

# mutant from group 1 
cp /home/admin4/2025_yeast_NGS/group_1_S1_L001_R1_001.fastq.gz .

# check quality of fastq files
fastqc /home/group_2/uz_real_data/1_raw_data/group_2_S2_L001_R1_001.fastq.gz
fastqc /home/group_2/uz_real_data/1_raw_data/control1_R1.fastq


# make index for the reference genome
hisat2-build /home/group_2/uz_real_data/1_raw_data/ncbi_dataset/data/GCF_009730115.1/GCF_009730115.1_ASM973011v1_genomic.fna sh

# align reads to the reference genome
hisat2 -q -x ./yeast -U control1_R1.fastq -S ../3_control_processing/control.sam > ../3_control_processing/align_control.log 2>&1
hisat2 -q -x ./yeast -U group_2_S2_L001_R1_001.fastq.gz -S ../3_mutant_processing_group2/mutant_group2_yeast.sam > ../3_mutant_processing_group2/align_group2_yeast.log 2>&1
hisat2 -q -x ./Staphylococcus_hominis/sh -U group_2_S2_L001_R1_001.fastq.gz -S ../3_mutant_processing_group2/mutant_group2_sh.sam > ../3_mutant_processing_group2/align_group2_sh.log 2>&1
hisat2 -q -x ./yeast -U group_1_S1_L001_R1_001.fastq.gz -S ../3_mutant_processing_group1/mutant_group1.sam > ../3_mutant_processing_group1/align_group1.log 2>&1

uz_real_data/1_raw_data/ Staphylococcus_hominis

# check the first 30 lines of the output
head -30 ym-WT.sam

# convert to bam
samtools view -bS -h 3_control_processing/control.sam > 3_control_processing/control.bam
samtools view -bS -h 3_mutant_processing_group1/mutant_group1.sam > 3_mutant_processing_group1/mutant_group1.bam

# visualize bam file
samtools view 3_mutant_processing_group1/mutant_group1.bam

# filter out unmapped reads
samtools view -b -F 4 3_control_processing/control.bam > 3_control_processing/control_mapped.bam
samtools view -b -F 4 3_mutant_processing_group1/mutant_group1.bam > 3_mutant_processing_group1/mutant_group1_mapped.bam

# sort bam file
gatk SortSam -INPUT 3_control_processing/control_mapped.bam -OUTPUT 3_control_processing/control_sorted.bam -SORT_ORDER coordinate
gatk SortSam -INPUT 3_mutant_processing_group1/mutant_group1_mapped.bam -OUTPUT 3_mutant_processing_group1/mutant_group1_sorted.bam -SORT_ORDER coordinate

# mark duplicates
gatk MarkDuplicates -INPUT 3_control_processing/control_sorted.bam -OUTPUT 3_control_processing/control_sorted_and_marked_file.bam -METRICS_FILE 3_control_processing/control_sorted_and_marked_file_metrics.txt -VALIDATION_STRINGENCY LENIENT -REMOVE_DUPLICATES True -CREATE_INDEX True
gatk MarkDuplicates -INPUT 3_mutant_processing_group1/mutant_group1_sorted.bam -OUTPUT 3_mutant_processing_group1/mutant_group1_sorted_and_marked_file.bam -METRICS_FILE 3_mutant_processing_group1/mutant_group1_sorted_and_marked_file.txt -VALIDATION_STRINGENCY LENIENT -REMOVE_DUPLICATES True -CREATE_INDEX True

# add read groups
gatk AddOrReplaceReadGroups -I 3_control_processing/control_sorted_and_marked_file.bam -O 3_control_processing/control_sorted_marked_ReadGroups_file.bam -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20
gatk AddOrReplaceReadGroups -I 3_mutant_processing_group1/mutant_group1_sorted_and_marked_file.bam -O 3_mutant_processing_group1/mutant_group1_sorted_marked_ReadGroups_file.bam -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20

# index bam file
samtools index 3_control_processing/control_sorted_marked_ReadGroups_file.bam
samtools index 3_mutant_processing_group1/mutant_group1_sorted_marked_ReadGroups_file.bam

# move to the folder with yeast genome
cd /home/group_2/uz_real_data/1_raw_data

# create index file for the reference genome with picard
PICARD=/home/admin4/SNP_calling/first_part/software/picard.jar
java -jar $PICARD CreateSequenceDictionary R=S288C_reference_plus_2u.fa O=S288C_reference_plus_2u.dict
samtools faidx S288C_reference_plus_2u.fa

cd ..

# Recalibrating base scores
gatk HaplotypeCaller -R 1_raw_data/S288C_reference_plus_2u.fa -I 3_control_processing/control_sorted_marked_ReadGroups_file.bam -O 3_control_processing/control_raw_variants.vcf > logs/control_raw_variants.log 2>&1
gatk SelectVariants -R 1_raw_data/S288C_reference_plus_2u.fa -V 3_control_processing/control_raw_variants.vcf --select-type-to-include SNP -O 3_control_processing/control_raw_snps.vcf > logs/control_raw_snps.log 2>&1
gatk SelectVariants -R 1_raw_data/S288C_reference_plus_2u.fa -V 3_control_processing/control_raw_variants.vcf --select-type-to-include INDEL -O 3_control_processing/control_raw_indels.vcf > logs/control_raw_indels.log 2>&1
gatk VariantFiltration -R 1_raw_data/S288C_reference_plus_2u.fa -V 3_control_processing/control_raw_snps.vcf --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filter-name "basic_snp_filter" -O 3_control_processing/control_filtered_snps.vcf > logs/control_filtered_snps.log 2>&1
gatk VariantFiltration -R 1_raw_data/S288C_reference_plus_2u.fa -V 3_control_processing/control_raw_indels.vcf --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filter-name "basic_indel_filter" -O 3_control_processing/control_filtered_indels.vcf > logs/control_filtered_indels.log 2>&1
gatk BaseRecalibrator -R 1_raw_data/S288C_reference_plus_2u.fa -I 3_control_processing/control_sorted_marked_ReadGroups_file.bam  --known-sites 3_control_processing/control_filtered_snps.vcf  --known-sites 3_control_processing/control_filtered_indels.vcf -O 3_control_processing/control_recal_data.table > logs/control_recal_data.log 2>&1
gatk ApplyBQSR -R 1_raw_data/S288C_reference_plus_2u.fa -I 3_control_processing/control_sorted_marked_ReadGroups_file.bam --bqsr-recal-file 3_control_processing/control_recal_data.table -O 3_control_processing/control_recal_reads.bam > logs/control_recal_reads.log 2>&1

# recalibrate for the other sample as well
gatk HaplotypeCaller -R 1_raw_data/S288C_reference_plus_2u.fa -I 3_mutant_processing_group1/mutant_group1_sorted_marked_ReadGroups_file.bam -O 3_mutant_processing_group1/mutant_group1_raw_variants.vcf > logs/mutant_group1_raw_variants.log 2>&1
gatk SelectVariants -R 1_raw_data/S288C_reference_plus_2u.fa -V 3_mutant_processing_group1/mutant_group1_raw_variants.vcf --select-type-to-include SNP -O 3_mutant_processing_group1/mutant_group1_raw_snps.vcf > logs/mutant_group1_raw_snps.log 2>&1
gatk SelectVariants -R 1_raw_data/S288C_reference_plus_2u.fa -V 3_mutant_processing_group1/mutant_group1_raw_variants.vcf --select-type-to-include INDEL -O 3_mutant_processing_group1/mutant_group1_raw_indels.vcf > logs/mutant_group1_raw_indels.log 2>&1
gatk VariantFiltration -R 1_raw_data/S288C_reference_plus_2u.fa -V 3_mutant_processing_group1/mutant_group1_raw_snps.vcf --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filter-name "basic_snp_filter" -O 3_mutant_processing_group1/mutant_group1_filtered_snps.vcf > logs/mutant_group1_filtered_snps.log 2>&1
gatk VariantFiltration -R 1_raw_data/S288C_reference_plus_2u.fa -V 3_mutant_processing_group1/mutant_group1_raw_indels.vcf --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filter-name "basic_indel_filter" -O 3_mutant_processing_group1/mutant_group1_filtered_indels.vcf > logs/mutant_group1_filtered_indels.log 2>&1

gatk BaseRecalibrator -R 1_raw_data/S288C_reference_plus_2u.fa -I 3_mutant_processing_group1/mutant_group1_sorted_marked_ReadGroups_file.bam  --known-sites 3_mutant_processing_group1/mutant_group1_filtered_snps.vcf  --known-sites 3_mutant_processing_group1/mutant_group1_filtered_indels.vcf -O 3_mutant_processing_group1/mutant_group1_recal_data.table > logs/mutant_group1_recal_data.log 2>&1
gatk ApplyBQSR -R 1_raw_data/S288C_reference_plus_2u.fa -I 3_mutant_processing_group1/mutant_group1_sorted_marked_ReadGroups_file.bam --bqsr-recal-file 3_mutant_processing_group1/mutant_group1_recal_data.table -O 3_mutant_processing_group1/mutant_group1_recal_reads.bam > logs/mutant_group1_recal_reads.log 2>&1

# call variants
gatk HaplotypeCaller -R 1_raw_data/S288C_reference_plus_2u.fa -I 3_control_processing/control_recal_reads.bam  -O 3_control_processing/control_variants_recal.vcf > logs/control_variants_recal.log 2>&1
gatk HaplotypeCaller -R 1_raw_data/S288C_reference_plus_2u.fa -I 3_mutant_processing_group1/mutant_group1_recal_reads.bam  -O 3_mutant_processing_group1/mutant_group1_variants_recal.vcf > logs/mutant_group1_variants_recal.log 2>&1

# Check how many SNP and indels were called
head 3_mutant_processing_group1/mutant_group1_variants_recal.vcf

# instead of counting manually
grep -c ^chr 3_control_processing/control_variants_recal.vcf # 204
grep -c ^chr 3_mutant_processing_group1/mutant_group1_variants_recal.vcf    # 714

# check the number of low quality variants
grep -c LowQual 3_control_processing/control_variants_recal.vcf   # 1
grep -c LowQual 3_mutant_processing_group1/mutant_group1_variants_recal.vcf # 1

# get the actual variant
grep LowQual 3_control_processing/control_variants_recal.vcf 
grep LowQual 3_mutant_processing_group1/mutant_group1_variants_recal.vcf

# compare the files to report the positions which are common and unique in mutant and wildtype.
vcftools --vcf 3_mutant_processing_group1/mutant_group1_variants_recal.vcf --diff 3_control_processing/control_variants_recal.vcf  --diff-site --out 4_final/compared > 4_final/compared.log 2>&1

# extract mutations, specific to mutant (printed in variants.txt)
less 4_final/compared.diff.sites_in_files | awk 'BEGIN{print "Chr\tPos\tREF\tALT"}; NR > 1{ if ($4 == 1) { print $1, $2, $5, $7} }' OFS='\t' > 4_final/variants.txt

