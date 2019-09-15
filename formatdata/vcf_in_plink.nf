#!/usr/bin/env nextflow

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
allowed_params = ['file_listvcf', 'min_scoreinfo']


params.help = false

nextflowversion =nextflow.version


if (workflow.repository)
  wflowversion="${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
else
  wflowversion="A local copy of the workflow was used"

report = new LinkedHashMap()

// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
}

params.queue    = 'batch'
params.remove_on_bp  = 1

params.file_listvcf=""
params.min_scoreinfo=0.6
params.max_plink_cores = 8 
params.plink_mem_req = '10GB' // how much plink needs for this
params.output_pat="out"
params.output_dir="plink/"



if(params.file_listvcf==""){
error('params.file_listvcf : file contains list vcf not found')
}
/*read file to have list of channel for each vcf*/
list_vcf=Channel.fromPath(file(params.file_listvcf).readLines())

process formatvcf{
  cpus params.max_plink_cores
  memory params.plink_mem_req
  time   params.big_time
  input :
     file(vcf) from list_vcf
  output :
     set file("${Ent}.bed"), file("${Ent}.bim"), file("${Ent}.fam") into listchroplink
     file("${Ent}.bim") into listbimplink
  script :
    Ent=vcf.baseName
    """
    bcftools view -Ou -i 'INFO>${params.min_scoreinfo}' $vcf  | bcftools convert -Oz -o ${Ent}.vcf.gz
    plink --vcf ${Ent}.vcf.gz --recode --keep-allele-order --make-bed --out ${Ent} --threads ${params.max_plink_cores}
    cp ${Ent}.bim ${Ent}.save.bim
    awk \'{if(\$2==\".\"){\$2=\$1\":\"\$4};print \$0}\' ${Ent}.save.bim > ${Ent}.bim
    """
}
bimmerg=listbimplink.collect()
process GetRsDup{
    input :
      file(bim) from bimmerg
    output :
      file(out) into duplicat
    script :
      lbim=bim.join(",")
      out="snpfile_red.rs"
      """
      search_dup_bim.py $lbim $out
      """ 
}
listchroplinkinf=duplicat.combine(listchroplink)
process TransformRsDup{
    input :
     set file(rstochange), file(bed),file(bim),file(fam) from listchroplinkinf
    output :
       set file(bed),file(bim),file(fam) into listchroplinkrs
       val(header) into plinkhead
    script :
       header=bed.baseName 
       """
       cp $bim ${bim}.save
       replacers_forbim.py ${bim}.save $rstochange $bim  
       """
}
listplinkinf=listchroplinkrs.collect()
headplink=plinkhead.collect()

process MergePlink{
  cpus params.max_plink_cores
  memory params.plink_mem_req
  time   params.big_time
  input :
       file(lplk) from listplinkinf
       val(hplk) from headplink
  publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
  output :
     set file("${params.output_pat}.bed"), file("${params.output_pat}.bim"),file("${params.output_pat}.fam") into plinkformatf
  script :
       hplkFirst=hplk[0]
       hplk.remove(0)
       hplk2=hplk.join(',')
       """
       echo $hplk2|awk -F\',\' \'{for(Cmt=1;Cmt<=NF;Cmt++)print \$Cmt"\\n"}\' > fileplk
       plink --bfile $hplkFirst --keep-allele-order --threads ${params.max_plink_cores} --merge-list fileplk --make-bed --out ${params.output_pat}
       """
}



