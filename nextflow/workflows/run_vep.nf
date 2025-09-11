/* 
 * Workflow to run VEP on multiple input files
 *
 * Requires Nextflow: https://nextflow.io
 */

nextflow.enable.dsl=2

// module imports
include { checkVCFheader; checkVCF } from '../modules/check_VCF.nf'
include { generateSplits } from '../modules/generate_splits.nf'
include { splitVCF } from '../modules/split_VCF.nf' 
include { mergeVCF } from '../modules/merge_VCF.nf'  
include { runVEP as runVEPonVCF } from '../modules/run_vep.nf'
include { runVEP } from '../modules/run_vep.nf'
include { filter_vep } from '../modules/filter_vep.nf'
include { split_snvs_indels } from '../modules/split_snvs_indels.nf'

// print usage
if (params.help) {
  log.info """
Pipeline to run VEP
-------------------

Usage:
  nextflow run main.nf -profile <profiles-to-use> --input <path-to-file> --vep_config <path-to-vep-config>

Parameters:

[Mandatory]
  --input FILE              Input file (if unsorted, use --sort to avoid errors in indexing the output file). Alternatively, can also be a directory containing input files
  --vep_config FILENAME     VEP config file. Alternatively, can also be a directory containing VEP INI files. Default: 'vep_config/vep.ini'

[OPTIONAL]
  --bin_size INT            Number of lines to split input into multiple jobs. Default: 100
  --vep_version VERSION     VEP version to use from Docker Hub (such as 113.0); only required when using Docker or Singularity profile. Default: 'latest'
  --cpus INT                Number of CPUs to use. Default: 1
  --outdir DIRNAME          Name of output directory. Default: outdir
  --output_prefix PREFIX    Output filename prefix. The generated output file will have name <output_prefix>_VEP.vcf.gz.
                            NOTE: Do not use this parameter if you are expecting multiple output files.

  --sort                    Sort VCF results from VEP (only required if input is unsorted; slower if enabled). Default: false
  --filters STRING          Comma-separated list of filter conditions to pass to filter_vep,
                            such as "AF < 0.01,Feature is ENST00000377918".
                            Read more on how to write filters at https://ensembl.org/info/docs/tools/vep/script/vep_filter.html
                            Default: null (filter_vep is not run)
  """
  exit 1
}

workflow run_vep {
  take:
    inputs
  main:


    runVEPonVCF(inputs)
    //// Run VEP on VCF files with header
    //inputs |
    //  checkVCF |
    //  // Generate split files that each contain bin_size number of variants
    //  generateSplits | transpose |
    //  // Split VCF using split files
    //  splitVCF | transpose |
    //  // Run VEP for each split VCF file and for each VEP config
    //  map { it + [format: 'vcf'] } | runVEPonVCF

    //// Merge split VCF files (creates one output VCF for each input VCF)
    //// COULD OPTIMIZE HERE (currently waits for all of above to finish)
    //out = runVEPonVCF.out.files
    //        .groupTuple(by: [0, 1, 4])
    //mergeVCF(out)
    filter_vep(runVEPonVCF.out)

    split_snvs_indels(filter_vep.out.vcf)

  emit:
    split_snvs_indels.out.snvs
}
