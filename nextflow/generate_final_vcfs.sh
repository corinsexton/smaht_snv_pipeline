#!/bin/bash

# Usage: $0 -o output.vcf.gz input1.vcf.gz [input2.vcf.gz ...]
./final_sets.sh -o ST002-1G_4callers.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/1_MT_to_ST002-1G_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/10_RUFUS_ST002-1G_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/5_STRELKA_ST002-1G_ST001-1D_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/7_longcallD_ST002-1G.phased.vcf.gz

./final_sets.sh -o ST002-1G_3callers.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/1_MT_to_ST002-1G_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/10_RUFUS_ST002-1G_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/5_STRELKA_ST002-1G_ST001-1D_BCM_300x.phased.vcf.gz 

./final_sets.sh -o ST002-1G_7callers.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/1_MT_to_ST002-1G_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/10_RUFUS_ST002-1G_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/5_STRELKA_ST002-1G_ST001-1D_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/2_MT_pu_ST002-1G_ST001-1D_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/9_deepsomatic_pmm_ST002-1G.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/6_deepsomatic_ST002-1G_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/7_longcallD_ST002-1G.phased.vcf.gz


./final_sets.sh -o ST002-1D_3callers.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/1_MT_to_ST002-1D_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/10_RUFUS_ST002-1D_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/2_STRELKA_ST002-1D_ST001-1D_BCM_300x.phased.vcf.gz 

#./final_sets.sh -o ST002-1D_4callers.vcf.gz \
#	pooled_truthset_PTA_results/phasing/step5/passed/1_MT_to_ST002-1D_BCM_300x.phased.vcf.gz \
#	pooled_truthset_PTA_results/phasing/step5/passed/8_RUFUS_ST002-1D_ST001-1D_BCM_300x.phased.vcf.gz \
#	pooled_truthset_PTA_results/phasing/step5/passed/2_STRELKA_ST002-1D_ST001-1D_BCM_300x.phased.vcf.gz \
#	pooled_truthset_PTA_results/phasing/step5/passed/7_longcallD_ST002-1D.phased.vcf.gz

./final_sets.sh -o ST002-1D_4callers.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/1_MT_to_ST002-1D_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/10_RUFUS_ST002-1D_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/2_STRELKA_ST002-1D_ST001-1D_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/7_longcallD_ST002-1D.phased.vcf.gz

./final_sets.sh -o ST002-1D_7callers.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/1_MT_to_ST002-1D_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/10_RUFUS_ST002-1D_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/2_STRELKA_ST002-1D_ST001-1D_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/2_MT_pu_ST002-1D_ST001-1D_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/9_deepsomatic_pmm_ST002-1D.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/6_deepsomatic_ST002-1D_BCM_300x.phased.vcf.gz \
	pooled_truthset_PTA_results/phasing/step5/passed/7_longcallD_ST002-1D.phased.vcf.gz


# 2 callers:
./final_sets.sh -o final_vcfs/single_gcc/2_caller_test/SMHT005-3C_2gccs.vcf.gz \
	production_results/phasing/step5/passed/MT2_to-SMHT005-3C.phased.vcf.gz \
	production_results/phasing/step5/passed/strelka-SMHT005-3C.phased.vcf.gz

./final_sets.sh -o final_vcfs/single_gcc/2_caller_test/SMHT004-3C_2gccs.vcf.gz \
	production_results/phasing/step5/passed/MT2_to-SMHT004-3C.phased.vcf.gz \
	production_results/phasing/step5/passed/strelka-SMHT004-3C.phased.vcf.gz

./final_sets.sh -o final_vcfs/single_gcc/2_caller_test/SMHT007-3C_2gccs.vcf.gz \
	production_results/phasing/step5/passed/MT2_to-SMHT007-3C.phased.vcf.gz \
	production_results/phasing/step5/passed/strelka-SMHT007-3C.phased.vcf.gz

./final_sets.sh -o final_vcfs/single_gcc/2_caller_test/SMHT008-3C_2gccs.vcf.gz \
	production_results/phasing/step5/passed/MT2_to-SMHT008-3C.phased.vcf.gz \
	production_results/phasing/step5/passed/strelka-SMHT008-3C.phased.vcf.gz

./final_sets.sh -o final_vcfs/single_gcc/2_caller_test/SMHT009-3C_2gccs.vcf.gz \
	production_results/phasing/step5/passed/MT2_to-SMHT009-3C.phased.vcf.gz \
	production_results/phasing/step5/passed/strelka-SMHT009-3C.phased.vcf.gz

# single GCC exercises:
./final_sets.sh -o final_vcfs/single_gcc/SMHT005-3C_Broad.vcf.gz \
	single_gccs_results/phasing/step5/passed/MTto_SMHT005-3C_Broad.phased.vcf.gz \
	single_gccs_results/phasing/step5/passed/strelka_SMHT005-3C_Broad.phased.vcf.gz 

./final_sets.sh -o final_vcfs/single_gcc/SMHT005-3C_UW.vcf.gz \
	single_gccs_results/phasing/step5/passed/MTto_SMHT005-3C_UW.phased.vcf.gz \
	single_gccs_results/phasing/step5/passed/strelka_SMHT005-3C_UW.phased.vcf.gz 


./final_sets.sh -o final_vcfs/single_gcc/SMHT004-3C_BCM.vcf.gz \
	single_gccs_results/phasing/step5/passed/MTto_SMHT004-3C_BCM.phased.vcf.gz \
	single_gccs_results/phasing/step5/passed/strelka_SMHT004-3C_BCM.phased.vcf.gz 

./final_sets.sh -o final_vcfs/single_gcc/SMHT004-3C_UW.vcf.gz \
	single_gccs_results/phasing/step5/passed/MTto_SMHT004-3C_UW.phased.vcf.gz \
	single_gccs_results/phasing/step5/passed/strelka_SMHT004-3C_UW.phased.vcf.gz 


./final_sets.sh -o final_vcfs/single_gcc/SMHT007-3C_Broad.vcf.gz \
	single_gccs_results/phasing/step5/passed/MTto_SMHT007-3C_Broad.phased.vcf.gz \
	single_gccs_results/phasing/step5/passed/strelka_SMHT007-3C_Broad.phased.vcf.gz 

./final_sets.sh -o final_vcfs/single_gcc/SMHT007-3C_NYGC.vcf.gz \
	single_gccs_results/phasing/step5/passed/MTto_SMHT007-3C_NYGC.phased.vcf.gz \
	single_gccs_results/phasing/step5/passed/strelka_SMHT007-3C_NYGC.phased.vcf.gz 


./final_sets.sh -o final_vcfs/single_gcc/SMHT008-3C_WashU.vcf.gz \
	single_gccs_results/phasing/step5/passed/MTto_SMHT008-3C_WashU.phased.vcf.gz \
	single_gccs_results/phasing/step5/passed/strelka_SMHT008-3C_WashU.phased.vcf.gz 

./final_sets.sh -o final_vcfs/single_gcc/SMHT008-3C_BCM.vcf.gz \
	single_gccs_results/phasing/step5/passed/MTto_SMHT008-3C_BCM.phased.vcf.gz \
	single_gccs_results/phasing/step5/passed/strelka_SMHT008-3C_BCM.phased.vcf.gz 


./final_sets.sh -o final_vcfs/single_gcc/SMHT009-3C_WashU.vcf.gz \
	single_gccs_results/phasing/step5/passed/MTto_SMHT009-3C_WashU.phased.vcf.gz \
	single_gccs_results/phasing/step5/passed/strelka_SMHT009-3C_WashU.phased.vcf.gz 

./final_sets.sh -o final_vcfs/single_gcc/SMHT009-3C_BCM.vcf.gz \
	single_gccs_results/phasing/step5/passed/MTto_SMHT009-3C_BCM.phased.vcf.gz \
	single_gccs_results/phasing/step5/passed/strelka_SMHT009-3C_BCM.phased.vcf.gz 
