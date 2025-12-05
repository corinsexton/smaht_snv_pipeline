#!/bin/bash

./minipileup-parallel_sr_only.sh \
	-i ../final_vcfs/3_callers/SMHT008/SMHT008-3AD.vcf.gz \
	-r ~/hg38_no_alt.fa -o output_test \
	--sr-cram /home/cos689/smaht/DATA/GCC_BCM/SMHT008/SMHT008-3AD/illuminaNovaseq_bulkWgs/seq_data/SMAFI8T6MYKW.bam \
	--sr-tissue 3AD \
	--sr-cram /home/cos689/smaht/DATA/GCC_BCM/SMHT008/SMHT008-3AF/illuminaNovaseq_bulkWgs/seq_data/SMAFILEXSW5I.bam \
	--sr-tissue 3AF -t 2
