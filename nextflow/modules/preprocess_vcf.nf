
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
	bcftools norm --check-ref x -m- -f ${ref} ${vcf} -Ou \
	  | bcftools view -f PASS -Oz -o ${id}.norm.PASS.vcf.gz -Wtbi
	bcftools norm -a -Oz -o ${id}.norm.PASS.atom.vcf.gz ${id}.norm.PASS.vcf.gz -Wtbi
	

    py_dedup_atomized.py ${id}.norm.PASS.vcf.gz ${id}.norm.PASS.atom.vcf.gz ${id}.norm.PASS.atom.dedup.vcf.gz
    """
}

