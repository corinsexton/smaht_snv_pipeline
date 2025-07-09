#!/bin/bash


#   ls /n/data1/hms/dbmi/park/jiny/SMaHT/Tissue6/truthset/1.LR_valid/*/output/*OUT* > all_donors_all_software_vcfs.txt
#   grep ST004-1Q all_donors_all_software_vcfs.txt | grep validOUT > ST004-1Q_vcfs.tsv
#   grep ST003-1Q all_donors_all_software_vcfs.txt | grep validOUT > ST003-1Q_vcfs.tsv
#   grep ST002-1G all_donors_all_software_vcfs.txt | grep validOUT > ST002-1G_vcfs.tsv
#   grep ST002-1D all_donors_all_software_vcfs.txt | grep validOUT > ST002-1D_vcfs.tsv
#   grep ST001-1D all_donors_all_software_vcfs.txt | grep validOUT > ST001-1D_vcfs.tsv
#   grep ST001-1A all_donors_all_software_vcfs.txt | grep validOUT > ST001-1A_vcfs.tsv


#for i in ST*tsv; do
#	id=${i%_vcfs.txt}
#	echo ${id}
#	while read -r f ; do
#		f=${f%.gz}
#		f_local=${f##*/}
#		./fix_malformed_vcfs.py -i ${f} -o bgzipped/prelim/${f_local}
#		bcftools view -G -Ob bgzipped/prelim/${f_local} -Wtbi -o bgzipped/${f_local}.gz
#	done < ${i}
#	echo finished bgzipping
#done

for tissue in ST004-1Q ST003-1Q ST002-1G ST002-1D ST001-1D ST001-1A; do
	ls bgzipped/*${tissue}*vcf.gz > ${tissue}_vcfs.tsv
done

for tissue in ST002 ST001; do
	ls bgzipped/*${tissue}*vcf.gz > ${tissue}_vcfs.tsv
done

for tissue in ST; do
	ls bgzipped/*${tissue}*vcf.gz > ${tissue}_vcfs.tsv
done

for f in *_vcfs.tsv; do
	echo Starting ${f}
	id=${f%_vcfs.tsv}
	echo ${id}
	bcftools merge --force-samples -m none \
		 -Ob -o concat_vcfs/panel_of_errors.${id}.vcf.gz -Wtbi \
		--file-list ${f}
done

cp concat_vcfs/panel_of_errors.ST.vcf.gz POE_benchmarking.vcf.gz
cp concat_vcfs/panel_of_errors.ST.vcf.gz.tbi POE_benchmarking.vcf.gz.tbi
