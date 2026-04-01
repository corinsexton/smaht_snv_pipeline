/* 
 * Workflow to run VEP on multiple input files
 *
 * Requires Nextflow: https://nextflow.io
 */

nextflow.enable.dsl=2

// module imports
include { runVEP } from '../modules/run_vep.nf'
include { filter_vep } from '../modules/filter_vep.nf'

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
    regions_input

  main:


    runVEP(inputs,'mosaic')
    filter_vep(runVEP.out, 'mosaic')

  emit:
    filter_vep.out.vcf
}
