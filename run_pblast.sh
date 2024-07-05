#!/usr/bin/bash

WD=/gpfs/gsfs12/users/lorenziha/KAREN_USDIN/TK_119/USDIN_LAB_TK_32/simulated_genome
PROBES=${WD}/exome_probes.mm9.fasta
OUTPUT=${WD}/sim_diployd_genome.psl
GENOME2BIT=${WD}/sim_diployd_genome.fa.2bit

source myconda

conda activate read_simulation 

pblat ${GENOME2BIT} ${PROBES} ${OUTPUT} -minScore=95 -minIdentity=95 -threads=8
