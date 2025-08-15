
// norm, filter by PASS, atomize, and remove any duplicates from SNV caller vcf files.

process preprocess_vcf {
    cpus 1
    memory '2G'
    time '30m'

    input:
    tuple val(id), path(vcf), path(tbi), path(ref)

    output:
    tuple val(id), path("${id}.norm.PASS.atom.dedup.vcf.gz"), path("${id}.norm.PASS.atom.dedup.vcf.gz.tbi")

    script:
    """
    bcftools view --header-only ${vcf} > header.txt

    if grep -q -m1 'FEX' header.txt; then
        # this is only for RUFUS files (FEX=PASS)
	    bcftools norm --check-ref x -m- -f ${ref} ${vcf} -Ou \
          | bcftools view -i 'FILTER=="PASS" || INFO/FEX == "PASS"' -Oz -o ${id}.norm.PASS.vcf.gz -Wtbi
	    bcftools norm -a -Oz -o ${id}.norm.PASS.atom.vcf.gz ${id}.norm.PASS.vcf.gz -Wtbi
    else
	    bcftools norm --check-ref x -m- -f ${ref} ${vcf} -Ou \
          | bcftools view -i 'FILTER=="PASS"' -Oz -o ${id}.norm.PASS.vcf.gz -Wtbi
	    bcftools norm -a -Oz -o ${id}.norm.PASS.atom.vcf.gz ${id}.norm.PASS.vcf.gz -Wtbi
    fi
	

    py_dedup_atomized.py ${id}.norm.PASS.vcf.gz ${id}.norm.PASS.atom.vcf.gz ${id}.norm.PASS.atom.dedup.vcf.gz
    """
}

