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


// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}



def helps = [ 'help' : 'help' ]
allowed_params = ['input_dir', 'input_pat', 'file_ref_gzip']

//params.file_vcf=""
params.input_dir="$PWD"
params.input_pat=""
params.file_ref_gzip=""
params.deleted_notref = "T"
params.output= "out"
params.output_dir= "out"

params.poshead_chro_inforef=0
params.parralchro = 0
params.poshead_bp_inforef=1
params.poshead_rs_inforef=2
params.poshead_a1_inforef=3
params.poshead_a2_inforef=4
params.bin_bcftools="bcftools"
params.bin_samtools="samtools"

params.plink_mem_req="10GB"
params.max_plink_cores="5"


inpat = "${params.input_dir}/${params.input_pat}"
println inpat
plkinit=Channel.create()
biminitial=Channel.create()
biminitial_extractref=Channel.create()
bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()

Channel
   .fromFilePairs("${inpat}.{bed,bim,fam}",size:3, flat : true){ file -> file.baseName }  \
      .ifEmpty { error "No matching plink files" }        \
      .map { a -> [checker(a[1]), checker(a[2]), checker(a[3])] }\
      .separate(plkinit, biminitial, biminitial_extractref) { a -> [a,a[1], a[1]] }


rs_infogz=Channel.fromPath(params.file_ref_gzip,checkIfExists:true)
process extractrsname{
  input :
    file(bim) from biminitial
    file(rsinfo) from rs_infogz
  output :
    file(outrs) into (l_infors1,l_infors2) 
  script :
    outrs="tmps_rsinfo"
    """
    awk '{print \$1" "\$4" "\$5" "\$6}' $bim > temp_pos
    zcat $rsinfo | extractrsid_bypos.py --file_chrbp temp_pos --out_file $outrs --ref_file stdin --chro_ps ${params.poshead_chro_inforef} --bp_ps ${params.poshead_bp_inforef} --rs_ps ${params.poshead_rs_inforef} --a1_ps ${params.poshead_a1_inforef}  --a2_ps ${params.poshead_a2_inforef}
    """
}

fastaextractref=Channel.fromPath(params.reffasta,checkIfExists:true)
process extractpositionfasta{
 input :
    file(bim) from biminitial_extractref
    file(fasta) from fastaextractref
  script :
    """
    extract_ref_bimf.py --bim $bim --fasta $fasta --out tmp
    """


}

process convertrsname{
  memory params.plink_mem_req
  cpus params.max_plink_cores
  input :
    file(rsinfo) from l_infors1
    set file(bed), file(bim), file(fam) from plkinit
  output :
    set file("${out}.bed"),file("${out}.bim"),file("${out}.fam") into plk_newrs
  script :
  tmpnewrs="newnamers"
  plk=bed.baseName
  out=plk+"_newrs"
  tmprs="filers"
  tmprsnodup="filers_nodup"
  extract=params.deleted_notref=='F' ? "" : " --extract range keep"
  """
  awk '{print \$1"\t"\$4}' $rsinfo > $tmprs
  biv_del_dup.py $tmprs  $tmprsnodup
  #awk '{print \$5}' $rsinfo > keep
  awk '{print \$1"\t"\$2"\t"\$2"\t"\$5}' $rsinfo > keep
  plink --keep-allele-order --noweb $extract --bfile $plk --make-bed --out $out --update-name $tmprsnodup -maf 0.0000000000000000001 --threads ${params.max_plink_cores}
  """
}

//remove multi allelic snps:

process deletedmultianddel{
   memory params.plink_mem_req
   cpus params.max_plink_cores
   input :
    set file(bed), file(bim), file(fam) from plk_newrs
   output :
    set file("${out}.bed"),file("${out}.bim"),file("${out}.fam") into plk_noindel
   script :
    plk=bed.baseName
    out=plk+"_nomulti"
    """
    biv_selgoodallele.r $bim rstodel
    plink --keep-allele-order --make-bed --bfile $plk --out $out -maf 0.0000000000000000001 --exclude rstodel --threads ${params.max_plink_cores}
    """  
}

process refallele{
   memory params.plink_mem_req
   cpus params.max_plink_cores
   input :
    set file(bed), file(bim), file(fam) from plk_noindel
    file(infors) from l_infors2
   output :
    set file("${out}.bed"),file("${out}.bim"),file("${out}.fam") into (plk_alleleref, plk_chrocount)
  script :
    plk=bed.baseName
    out=plk+"_refal"
    """
    awk '{print \$5"\t"\$6}' $infors > alleref
    plink --bfile $plk --make-bed --out $out --threads ${params.max_plink_cores} --a2-allele alleref
    """
}

hgrefconv=Channel.fromPath(params.reffasta,checkIfExists:true)
if(params.parralchro==0){
process convertInVcf {
   memory params.plink_mem_req
   cpus params.max_plink_cores
   time params.big_time
   input :
    set file(bed), file(bim), file(fam) from plk_alleleref
    file(fast) from hgrefconv
   publishDir "${params.output_dir}/vcf/", overwrite:true, mode:'copy'
   output :
    file("reportfixref")
    file("${out}.vcf.gz")  into (vcfi, vcfi2)
   script:
     base=bed.baseName
     out="${params.output}"
     """
     mkdir -p tmp
     plink --bfile ${base}  --recode vcf --out $out --keep-allele-order --snps-only --threads ${params.max_plink_cores} 
     ${params.bin_bcftools} sort  ${out}.vcf -O z > ${out}_tmp.vcf.gz -T tmp/
     ${params.bin_bcftools} +fixref ${out}_tmp.vcf.gz -Oz -o ${out}.vcf.gz -- -f $fast T tmp/ -m top &> reportfixref
     rm -rf tmp
     """
 }
}else{

process CounChro{
  input :
       set file(bed), file(bim), file(fam) from plk_chrocount
  output :
        stdout into chrolist
  script :
    """
    awk '{print \$1}' $bim|uniq 
    """
}

check2 = Channel.create()
ListeChro=chrolist.flatMap { list_str -> list_str.split() }.tap ( check2)

process convertInVcfChro{
   memory params.plink_mem_req
   cpus params.max_plink_cores
   input :
      set file(bed), file(bim), file(fam) from plk_alleleref
      file(fast) from hgrefconv
   publishDir "${params.output_dir}/vcf_bychro/", overwrite:true, mode:'copy'
   each chro from ListeChro
   output : 
     file("${out}.rep") 
     file("${out}.vcf.gz") into vcf_chro
   script:
     base=bed.baseName
     out="${params.output}"+"_"+chro
     """
     plink  --chr $chro --bfile ${base}  --recode vcf --out $out --keep-allele-order --snps-only --threads ${params.max_plink_cores}
     ${params.bin_bcftools} sort  ${out}.vcf -O z > ${out}_tmp.vcf.gz
     ${params.bin_bcftools} +fixref ${out}_tmp.vcf.gz -Oz -o ${out}.vcf.gz -- -f $fast -m top &> $out".rep"
     """
}

vcf=vcf_chro.collect()


process mergevcf{
  cpus params.max_plink_cores
  input :
   val(allfile) from vcf   
  publishDir "${params.output_dir}/vcf/", overwrite:true, mode:'copy'
  output :
     file("${out}.vcf.gz")  into (vcfi, vcfi2)
  script :
    fnames = allfile.join(" ")
    out="${params.output}"
    """  
    bcftools concat -Oz -o ${out}.vcf.gz --threads ${params.max_plink_cores} $fnames
    """
}

}

//if(params.reffasta!=""){
hgref=Channel.fromPath(params.reffasta, checkIfExists:true)
hgref2=Channel.fromPath(params.reffasta, checkIfExists:true)
process checkfixref{
  input :
    file(vcf) from vcfi
    file(hg) from hgref
  publishDir "${params.output_dir}/check/Bcftools", overwrite:true, mode:'copy'
  output :
    file("${params.output}.checkbcf*") 
  script :
    """
    ${params.bin_samtools} faidx $hg
    ${params.bin_bcftools} +fixref $vcf -- -f $hg 1> ${params.output}".checkbcf.out" 2> ${params.output}".checkbcf.err"
    """
}

process checkVCF{
  input :
    file(vcf) from vcfi2
    file(hg) from hgref2
  publishDir "${params.output_dir}/check/CheckVCF", overwrite:true, mode:'copy'
  output :
    file("${out}*") 
  script :
    out="${params.output}_check"
    """
    ${params.bin_samtools} faidx $hg
    checkVCF.py -r $hg -o $out $vcf
    """
}

//}
