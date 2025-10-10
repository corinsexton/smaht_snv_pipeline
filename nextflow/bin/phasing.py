#!/bin/bash


# step1: isolate only tier1 variants
# step2: find regions where to call haplotypecaller (7.5k around variants)
#   - py_variant_in_bed_7.5kflank.py
#   - merge intervals
# step3: run HaplotypeCaller
# step4: run vep again, keep anything > 0.001
#   - py_dbSNPf_common.py
#   - py_pos_for_phasing.py
#   - get_reads_for_phasing.py
#   - Filter_phasing4.py
