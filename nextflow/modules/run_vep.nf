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
  tuple val(id), val(original_file), path(index), path(vep_config)
  
  output:
  tuple val(id), path("${out}{.gz,}"), path("${out}{.gz,}.{tbi,csi}"), emit: files

  script:
  index_type = "tbi"
  out = file(original_file).getSimpleName() + "_VEP.vcf"
  tabix_arg = index_type == 'tbi' ? '' : '-C'
  
  vep_cmd = """
            vep --max_af -i ${original_file} \
            --custom file=/custom_ann/gnomad.joint.v4.1.sites.reduced.vcf.gz,short_name=gnomad4.1,format=vcf,type=exact,coords=0,fields=AF_grpmax_joint \
            --vcf --config ${vep_config}  -o out.vcf
            """
  

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
