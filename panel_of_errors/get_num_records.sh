#!/bin/bash

rm num_variants
for i in concat_vcfs/*tbi; do echo ${i} >> num_variants.tsv ; bcftools index --nrecords $i >> num_variants.tsv ; done
