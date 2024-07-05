# Generation of simulated genome and reads for evaluation of SNV calling pipelines  

### Concatenate SNP and INDEL data from SPRET_EiJ mouse strain
```
cat SPRET_EiJ.mgp.v5.indels.dbSNP142.normed.vcf.gz SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz > SPRET_EiJ.mgp.v5.ALL.vcf.gz
```

### Generate SNP and INDEL models
```
perl vcf2model.pl -vcf ./SPRET_EiJ.mgp.v5.ALL.vcf.gz -prefix SPRET_EiJ.SNP_INDEL_model
```

### Adding SNPs and INDELs anywhere
```
perl ~/data/bin/simuG/simuG.pl -refseq ../genome.fa -snp_count 10000 -indel_count 10000 -snp_model SPRET_EiJ.SNP_INDEL_model.SNP_model.txt -indel_model SPRET_EiJ.SNP_INDEL_model.INDEL_model.txt -prefix  mm10.10000snp.10000indel.SPRET_EiJ.model  -seed 1785035518
```

### Coding regions only
```
perl ~/data/bin/simuG/simuG.pl -refseq ../genome.fa -snp_count 10000 -indel_count 10000 -snp_model SPRET_EiJ.SNP_INDEL_model.SNP_model.txt -indel_model SPRET_EiJ.SNP_INDEL_model.INDEL_model.txt -prefix  mm10.10000snp.10000indel.SPRET_EiJ.model  -seed 1785035518 -coding_partition_for_snp_simulation coding -gene_gff Mus_musculus.GRCm38.102.chr.gff3.gz
```

### Generate sim genome for read generation

half homozygote mutations / half heterozygote mutations  
```
seqtk subseq mm10.10000snp.10000indel.SPRET_EiJ.model.simseq.DIPLOID.genome.fa seq.ids | sed 's/^>chr/>CHR/' > mm10.10000snp.10000indel.SPRET_EiJ.model.simseq.DIPLOID.genome_1_of_2.fa
seqtk subseq ../genome.fa seq2.ids | sed 's/^>chr/>CHR/' > mm10.10000snp.10000indel.SPRET_EiJ.model.simseq.DIPLOID.genome_2_of_2.fa

cat mm10.10000snp.10000indel.SPRET_EiJ.model.simseq.DIPLOID.genome.fa mm10.10000snp.10000indel.SPRET_EiJ.model.simseq.DIPLOID.genome_1_of_2.fa mm10.10000snp.10000indel.SPRET_EiJ.model.simseq.DIPLOID.genome_2_of_2.fa > sim_diployd_genome.fa

./faToTwoBit -long ./sim_diployd_genome.fa ./sim_diployd_genome.fa.2bit
```

### Example

pblat Homo_sapiens_assembly38.fasta.2bit ../w-Wessim/real_wes_reads_probes.fa blatoutput_Homo_sapiens_assembly38.psl -minScore=95 -minIdentity=95 -threads=8

```
sbatch --time=3-00:00:00 --mem=64g --cpus-per-task=8 run_pblast.sh
```

### ==============
### Using SimuSCoP
### ==============

### Calling SNPs and INDELs from original WES data 
```
./call_snps_with_gatk.sh genome.fa LB2BAV7596 ./01-trimming/LB2BAV7596_R1_trimmed_qc.fastq.gz ./01-trimming/LB2BAV7596_R2_trimmed_qc.fastq.gz
```

### Generating target bed file
```
grep '^chr' Mus_musculus.GRCm38.102.chr.CHRid.gff3| perl -ne '@x=split /\t/;$x[3]--; next if $x[2] ne CDS;     print "$x[0]\t$x[3]\t$x[4]\n"'| sort -u | sort -k1,1 -k2,2n | bedtools merge -i - > target.bed
```

### Generate SNP and INDEL profile
```
module load samtools

./simuscop/bin/seqToProfile -b ../LB2BAV7596.sorted.bam -t target.bed -v ../ALL_filtered_snps_indels.vcf -k 4 -B 100 -r ../genome.fa -s /usr/local/apps/samtools/1.19/bin/samtools -o LB2BAV7596.profile
```

### Make variation file
```
cat mm10.10000snp.10000indel.SPRET_EiJ.model.refseq2simseq.map.txt |perl ./build_variants_file.pl > variation_file.txt
```

### Estimate median fragment size from real sequencing data
```
module load deeptools
bamPEFragmentSize --bamfiles ../LB2BAV7596.sorted.bam 
```

### Generate simulated reads

option 1:
```
./simuscop/bin/simuReads config_LB2BAV7596_wes.txt
```

option 2:
```
sbatch ./generate_simulated_reads.sh
```

### Subsample simulated fastq files to get approx 40M reads per file
```
time seqtk sample -s100 control_1.fq 0.03 > control_0.03.R1.fastq
time seqtk sample -s100 control_2.fq 0.03 > control_0.03.R2.fastq
```



