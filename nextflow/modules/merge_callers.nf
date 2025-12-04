
// norm, filter by PASS, atomize, and remove any duplicates from SNV caller vcf files.
// then merge all callers into a single vcf

process merge_callers {
    cpus 1
    memory '16G'
    time '30m'

    publishDir "${params.results_dir}/merged_vcf",
    pattern: "${id}.merged.*.tsv",
    mode:'copy'

    publishDir "${params.results_dir}/merged_vcf",
    pattern: "${id}.merged.vcf.gz*",
    mode:'copy'

    input:
        tuple val(id), 
          val(callers),      // list of caller names
          path(vcfs),         // list of VCFs
          path(tbis),         // list of TBIs
          path(truth_vcfs, stageAs: "?/*"),   // list of truth VCFs
          path(truth_tbis, stageAs: "?/*")    // list of truth TBIs


    output:
        tuple val(id), path("${id}.merged.vcf.gz"), path("${id}.merged.vcf.gz.tbi"), path("truth.vcf.gz"), path("truth.vcf.gz.tbi"), emit:vcf
        path("${id}.merged.metrics.tsv"), emit: metrics
        path("${id}.merged.regions.tsv"), emit: regions

    /*
     callers structure example:
     [
       [ "SMHT004", "LR",  "LR_processed.vcf.gz",  "LR_processed.vcf.gz.tbi", truth_vcf, truth_tbi ],
       [ "SMHT004", "SR",  "SR_processed.vcf.gz",  "SR_processed.vcf.gz.tbi", truth_vcf, truth_tbi ],
       [ "SMHT004", "ONT", "ONT_processed.vcf.gz", "ONT_processed.vcf.gz.tbi", truth_vcf, truth_tbi ]
     ]
    */

    script:

    // Extract VCF paths from the list
    def caller_vcf_pairs = (0..<callers.size()).collect { idx ->
        "-i ${callers[idx]}:${vcfs[idx]}"
    }.join(" ")

    def truth_vcf = truth_vcfs[0]
    def truth_tbi = truth_tbis[0]

    """
    cp ${truth_vcf} truth.vcf.gz
    cp ${truth_tbi} truth.vcf.gz.tbi

    echo "Merging callers for ${id}"

    merge_callers.py ${caller_vcf_pairs} -s ${id} -o ${id}.merged.vcf.gz

    # --- metrics (standard schema) ---
    AFTER_VCF=${id}.merged.vcf.gz

    num_before=NA
    num_after=\$(bcftools view -H "\${AFTER_VCF}"  | wc -l | awk '{print \$1}')

    touch ${id}.merged.regions.tsv


    # Compute truth overlaps only if truth files are present
    if [[ -f "${truth_vcf}" ]]; then
      num_truth_before=NA
      num_truth_after=\$( bcftools isec -n=2 -w1 -c both "\${AFTER_VCF}"  "${truth_vcf}" 2>/dev/null | grep -v '^#' | wc -l | awk '{print \$1}')
    else
      num_truth_before=NA
      num_truth_after=NA
    fi

    {
      echo -e "id\tstep\tnum_before\tnum_truth_before\tnum_after\tnum_truth_after"
      echo -e "${id}\tmerged_callers\t\${num_before}\t\${num_truth_before}\t\${num_after}\t\${num_truth_after}"
    } > ${id}.merged.metrics.tsv
    """
}


