#!/usr/bin/env nextflow

/* 
 * Script to run VEP on split VCF files
 */

nextflow.enable.dsl=2

process runVEP {
  cpus 1
  memory '4G'
  time '2h'
  /*
  Run VEP on VCF files

  Returns
  -------
  Tuple of original VCF, split VCF file after running VEP, tabix index of that file, vep config file, a output dir, and the index type of VCF file
  */
  
  label 'vep'

  input:
  tuple val(id), val(original_file), path(input), path(index), path(vep_config), val(format)
  
  output:
  tuple val(id), val(original_file), path("${out}{.gz,}"), path("${out}{.gz,}.{tbi,csi}"), val(vep_config), emit: files

  script:
  index_type = "tbi"
  out = "vep" + "-" + file(original_file).getSimpleName() + "-" + vep_config.getSimpleName() + "-" + input.getName().replace(".gz", "")
  tabix_arg = index_type == 'tbi' ? '' : '-C'
  
  if( !input.exists() ) {
    exit 1, "Missing input: ${input}"

  }
  else if ( format == 'vcf' && !index.exists() ){
    exit 1, "VCF index file is not generated: ${index}"
  }
  else {
  def filters = "AF < 0.001,..."
  def filter_arg = ""
  for (filter in filters) {
    filter_arg = filter_arg + "-filter \"" + filter + "\" "
  }
  vep_cmd = """
            vep --max_af -i ${input} \
            --custom file=/custom_ann/gnomad.joint.v4.1.sites.reduced.vcf.gz,short_name=gnomad4.1,format=vcf,type=exact,coords=0,fields=AF_grpmax_joint \
            --vcf --config ${vep_config} | \
            filter_vep -o out.vcf --force_overwrite --only_matched ${filter_arg}
            """
  }

  if( params.sort ) {
    sort_cmd = "(head -1000 out.vcf | grep '^#'; grep -v '^#' out.vcf | sort -k1,1d -k2,2n) > ${out}"
  } else {
    sort_cmd = "mv out.vcf ${out}"
  }


  """
  ${vep_cmd}

  # Sort, bgzip and tabix VCF
  ${sort_cmd}

  bgzip ${out}
  tabix ${tabix_arg} ${out}.gz
  """
}
