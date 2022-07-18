#!/usr/bin/env nextflow

/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      jean-tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2022
 *
 *
 * Description  : Nextflow pipeline to transform vcf file in plink and other format
 *
 *(C) University of the Witwatersrand, Johannesburg, 2016-2019 on behalf of the H3ABioNet Consortium
 *This is licensed under the MIT Licence. See the "LICENSE" file for details
 */

//---- General definitions --------------------------------------------------//


import java.nio.file.Paths;
import sun.nio.fs.UnixPath;
import java.security.MessageDigest;
nextflow.enable.dsl = 1



def helps = [ 'help' : 'help' ]
allowed_params = ['file_listvcf', 'min_scoreinfo', "output_dir", "max_cores", "output", "bgen_bits", "mem_req", "genotype_field", "qctoolsv2_bin", "bcftools_bin", "score_imp", "bgen_type"]

params.mem_req = '10GB' // how much plink needs for this
params.output_dir="bgen/"
params.output="bgen"
params.file_listvcf=""
params.min_scoreinfo=0.6
params.max_cores = 8
params.genotype_field="GP"
params.qctoolsv2_bin="qctool"
params.bcftools_bin="bcftools"
params.score_imp="INFO"
params.bgen_type="bgen"
params.other_opt=""
params.bgen_bits=8


if(params.file_listvcf==""){
error('params.file_listvcf : file contains list vcf not found')
}
list_vcf=Channel.fromPath(file(params.file_listvcf).readLines())

process filter_vcf{
  label 'py3utils'
  cpus params.max_cores
  memory params.mem_req
  time   params.big_time
  input :
     file(vcf) from list_vcf
  output :
     set env(chro), file("${Ent}"), file("${Ent}.csi") into list_vcf_filt
  script :
    Ent=vcf.baseName+"_filter.vcf.gz"
    """
    ${params.bcftools_bin} view -i '${params.score_imp}>${params.min_scoreinfo}' $vcf -Oz --threads  ${params.max_cores} > $Ent
    ${params.bcftools_bin} index $Ent
    chro=`zcat $vcf|grep -v "#"|head -1|awk '{print \$1}'`

    """
}

process formatvcfinbgen{
  label 'py3utils'
  time   params.big_time
  memory params.mem_req
  cpus params.max_cores
  publishDir "${params.output_dir}/bgen_chro", overwrite:true, mode:'copy'
  input :
     tuple val(chro),file(vcf), file("${Ent}.csi") from list_vcf_filt
  output :
     file("${out}.bgen") into list_bgen
     file("${out}.sample") into list_bgen_sample
   script :
    out="${params.output}_${chro}"
    """
    ${params.qctoolsv2_bin} -g $vcf -vcf-genotype-field ${params.genotype_field} -ofiletype ${params.bgen_type} -og ${out}.bgen -filetype vcf -os ${out}.sample ${params.other_opt} -bgen-bits ${params.bgen_bits}
    """
}
/*
process formatvcfinbgen{
  label 'py3utils'
  memory params.mem_req
  time   params.big_time
  input :
     file(vcf) from list_bgen
  publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
  output :
     file("${Ent}.bgen")
     file("${Ent}.sample")
  script :
    Ent="${params.output}"
    """
    ${params.qctoolsv2_bin} -g $vcf -vcf-genotype-field ${params.genotype_field} -ofiletype ${params.bgen_type} -og ${Ent}.bgen -filetype vcf -os ${Ent}.sample ${params.other_opt} -bgen-bits ${params.bgen_bits}
    """
}

*/
lbgen=list_bgen.collect()
lsample=list_bgen_sample.collect()
process format_mergebgen{
  memory params.mem_req
  time   params.big_time
  label 'py3utils'
    input :
       path(lbgen) from lbgen
       path(lbgen) from lsample
    publishDir "${params.output_dir}/bgen_chro", overwrite:true, mode:'copy'
    output :
       file("${Ent}.bgen")
       file("${Ent}.sample")
    script :
    Ent="${params.output}"
    """
    filesample=`ls *.sample|head -1`
    ${params.qctoolsv2_bin} -g ${params.output}_#.bgen -ofiletype ${params.bgen_type} -og ${Ent}.bgen -filetype bgen -bgen-bits ${params.bgen_bits}  -s \$filesample -os ${Ent}.sample
    """
}

