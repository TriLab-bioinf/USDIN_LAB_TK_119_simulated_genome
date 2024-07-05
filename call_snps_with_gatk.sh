#!/usr/bin/bash

set -e

JAR=/opt/anaconda3/envs/tk_119//opt/gatk-3.8/GenomeAnalysisTK.jar
TMPDIR=TMP
BWA_DB=$1 
OUTPUT=ALL
PREFIX=$2 # sample prefix
#BAM_FILES="$@"
FASTQ1=$3
FASTQ2=$4
BAM=${PREFIX}.sorted.bam

module load bwa GATK samtools 


# ## Step 1
# # Mapping trimmed reads to reference genome with BWA
echo; echo Mapping reads to ${BWA_DB}; echo
if [[ ! -f ${BAM} ]]; then
bwa mem -R "@RG\tID:${PREFIX}\tPL:Illumina\tSM:${PREFIX}\tLB:${PREFIX}" -t 8 ${BWA_DB} ${FASTQ1} ${FASTQ2} \
        | samtools view -@ 8 -b -u -F 4 - | samtools sort -@ 8 -T TMP -O BAM -o ${BAM} --write-index  -
fi

BAM_FILES=${BAM}

for BAM in ${BAM_FILES}
do
    # LB2BAV7596_sorted_RG.filter.bam
    echo Processing... $BAM
     

    PREFIX_ARRAY=("${PREFIX_ARRAY[@]}" ${PREFIX})

    if [[ ! -f ${PREFIX}.raw.snps.indels.g.vcf ]]; then
            echo; echo "Running gatk HaplotypeCaller -ERC GVCF on ${BAM}"; echo
            gatk --java-options "-Xmx8G" HaplotypeCaller \
            -R ${BWA_DB} \
            -I ${BAM} \
            -ERC GVCF \
            --stand-call-conf 30.0  \
            -O ${PREFIX}.raw.snps.indels.g.vcf \
            -create-output-variant-index true \
            --sample-ploidy 2
    fi
done


for PREFIX in "${PREFIX_ARRAY[@]}"; do
    merge_param+=" --variant ${PREFIX}.raw.snps.indels.g.vcf"
done

# ## Step 4
# # Merge all desired g.vcf files from CombineGVCFs into a single VCF file of all strains
# # You will need to change the output of this name as well as making sure the input matches what you created in Step 3
# # This can run for a VERY long time. 96 Toxo genomes ran for 12 days straight before finishing
#cd $workingdir
if [[ ! -f ${OUTPUT}_raw_variants.vcf ]]; then
echo; echo "Running GenotypeGVCFs on *.g.vcf files" ; echo
gatk --java-options "-Xmx8G -Djava.io.tmpdir=${TMPDIR} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true'" GenotypeGVCFs \
    -R ${BWA_DB} \
    --max-alternate-alleles 2 \
    --stand-call-conf 30.0 \
    -O ${OUTPUT}_raw_variants.vcf \
    ${merge_param}
fi

## Step 5
##1. Extract the SNPs from the call set
##Action

##Run the following GATK command:
if [[ ! -f ${OUTPUT}_raw_snps_indels.vcf ]]; then
echo; echo "Running SelectVariants on ${OUTPUT}_raw_variants.vcf" ; echo
gatk --java-options "-Xmx8G -Djava.io.tmpdir=${TMPDIR}" SelectVariants \
    -R ${BWA_DB} \
    -V ${OUTPUT}_raw_variants.vcf \
    --select-type-to-include SNP \
    --select-type-to-include INDEL \
    -O ${OUTPUT}_raw_snps_indels.vcf
fi

## This creates a VCF file called raw_snps.vcf, containing just the SNPs from the original file of raw variants.

## 2. Apply the filter to the SNP call set
## Action

## Run the following GATK command:
if [[ ! -f ${OUTPUT}_filtered_snps_indels.vcf ]]; then
echo; echo "Running VariantFiltration on ${OUTPUT}_raw_snps.vcf" ; echo
gatk --java-options "-Xmx8G -Djava.io.tmpdir=${TMPDIR}" VariantFiltration \
    -R ${BWA_DB} \
    -V ${OUTPUT}_raw_snps_indels.vcf \
    --filter-name  "filter_QD" --filter-expression "QD < 2.0" \
    --filter-name  "filter_FS" --filter-expression "FS > 60.0" \
    --filter-name  "filter_MQ" --filter-expression "MQ < 40.0" \
    --filter-name  "filter_MQRankSum" --filter-expression "MQRankSum < -12.5" \
    --filter-name  "ftr_ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" \
    -O ${OUTPUT}_filtered_snps_indels.vcf
fi

