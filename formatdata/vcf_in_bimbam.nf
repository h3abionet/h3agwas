
/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      jean-tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2019
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


def helps = [ 'help' : 'help' ]
allowed_params = ['file_listvcf', 'min_scoreinfo', "output_dir", "max_plink_cores"]

params.plink_mem_req = '10GB' // how much plink needs for this
params.output_dir="impute2/"
params.file_listvcf=""
params.min_scoreinfo=0.6
params.max_plink_cores = 8
params.genotype_field="GP"
params.qctoolsv2_bin="qctool_v2"
params.bcftools_bin="bcftools"
params.score_imp="INFO"



if(params.file_listvcf==""){
error('params.file_listvcf : file contains list vcf not found')
}

/*read file to have list of channel for each vcf*/
list_vcf=Channel.fromPath(file(params.file_listvcf).readLines())

/*#/opt/exp_soft/bioinf/bin/qctool_v2.0.1 -g ${PATH_ASSOC}/chr${i}.vcf.gz -vcf-genotype-field GP -ofiletype bimbam_dosage -og ${PATH_ASSOC}/chr${i}.bimbam

#file check step
#zgrep -v '#' ${PATH_ASSOC}/chr${i}.vcf.gz | wc -l > ${PATH_ASSOC}/chr${i}.vcf.count
#grep '^[rs\|.:]' ${PATH_ASSOC}/chr${i}.bimbam | wc -l > ${PATH_ASSOC}/chr${i}.bimbam.count
cmp -s ${PATH_ASSOC}/chr${i}.vcf.count ${PATH_ASSOC}/chr${i}.bimbam.count && echo "for chr${i} vcf and bim files are equal" >> ${PATH_ASSOC}/check_vcf_to_BIMBAM.txt
*/

process formatvcfinbimbam{
  cpus params.max_plink_cores
  memory params.plink_mem_req
  time   params.big_time
  input :
     file(vcf) from list_vcf
  publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
  output :
     file("${Ent}.bimbam") 
  script :
    Ent=vcf.baseName
    """
    ${params.bcftools_bin} view -i '${params.score_imp}>${params.min_scoreinfo}' $vcf |${params.qctoolsv2_bin} -g - -vcf-genotype-field ${params.genotype_field} -ofiletype bimbam_dosage -og ${Ent}.bimbam -filetype vcf
    """
}




