
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
allowed_params = ['file_gwas', 'file_ref_gzip', "output_dir","output", "input_dir", "input_pat"]
allowed_header = ['head_pval', 'head_freq', 'head_bp', 'head_chr', 'head_rs', 'head_beta', 'head_se', 'head_A1', 'head_A2', 'sep']
allowed_headnewernew = ['headnew_pval', 'headnew_freq', 'headnew_bp', 'headnew_chr', 'headnew_rs', 'headnew_beta', 'headnew_se', 'headnew_A1', 'headnew_A2']
//chro_ps 0 --bp_ps 1 --rs_ps 
allowed_posfilref=['poshead_chro_inforef', 'poshead_bp_inforef','poshead_rs_inforef']
allowed_params+=allowed_header
allowed_params+=allowed_headnewernew


params.mem_req = '10GB' // how much plink needs for this
params.file_gwas=""
params.output_dir="out"
params.output="output"
params.input_dir=""
params.input_pat=""

params.head_pval = ""
params.head_freq = ""
params.head_bp = ""
params.head_chr = ""
params.head_rs = ""
params.head_beta=""
params.head_se=""
params.head_A1=""
params.head_A2=""
params.head_N=""
params.sep=""

params.headnew_pval = ""
params.headnew_freq = ""
params.headnew_bp = ""
params.headnew_chr = ""
params.headnew_rs = "rs"
params.headnew_beta=""
params.headnew_se=""
params.headnew_A1=""
params.headnew_A2=""
params.headnew_N=""

params.poshead_chro_inforef=0
params.poshead_bp_inforef=1
params.poshead_rs_inforef=2
params.poshead_a1_inforef=3
params.poshead_a2_inforef=4

if(params.file_gwas==""){
error('params.file_gwas: file contains gwas not found')
}

if(params.headnew_pval=="")params.headnew_pval=params.head_pval 
if(params.headnew_freq=="")params.headnew_freq=params.head_freq 
if(params.headnew_bp=="")params.headnew_bp=params.head_bp 
if(params.headnew_chr=="")params.headnew_chr=params.head_chr 
if(params.headnew_beta=="")params.headnew_beta=params.head_beta 
if(params.headnew_se=="")params.headnew_se=params.head_se 
if(params.headnew_A1=="")params.headnew_A1=params.head_A1 
if(params.headnew_A2=="")params.headnew_A2=params.head_A2 
//if(params.headnew_N=="")params.headnew_N=params.head_N 


gwas_chrolist = Channel.fromPath(params.file_gwas)
gwas_chrolist_ext = Channel.fromPath(params.file_gwas)
gwas_chrolist_ext = Channel.fromPath(params.file_gwas)
process getListeChro{
        input :
          file(gwas_res) from gwas_chrolist
        output :
          file("filechro") into chrolist
        script:
         sep=(params.sep!="") ?  " --sep ${params.sep}" : ""
         """
         extractlistchro.py --input_file $gwas_res --chro_header ${params.head_chr} $sep > filechro
        """
}


chrolist2=Channel.create()
chrolist.flatMap { list_str -> list_str.readLines()[0].split() }.set { chrolist2 }

process ExtractChroGWAS{
    memory params.mem_req 
    input :
      file(gwas) from gwas_chrolist_ext
    each chro from  chrolist2
    output :
      set val(chro), file(gwas_out) into gwas_format_chro
    script :
      gwas_out=gwas.baseName+"_"+chro+".gwas"
      sep=(params.sep!="") ?  "" : ""
      infofile="Chro:${params.head_chr}:${params.headnew_chr},Pos:${params.head_bp}:${params.headnew_bp},A2:${params.head_A2}:${params.headnew_A2},A1:${params.head_A1}:${params.headnew_A1},af:${params.head_freq}:${params.headnew_freq},Beta:${params.head_beta}:${params.headnew_beta},Se:${params.head_se}:${params.headnew_se},Pval:${params.head_pval}:${params.headnew_pval},N:${params.head_N}:${params.headnew_N}"
      """
      extractandformat_gwas.py --input_file $gwas --out_file ${gwas_out} --chr $chro --info_file $infofile
      """
}
gwas_format_chro_rs=gwas_format_chro.combine(Channel.fromPath(params.file_ref_gzip))

process ExtractRsIDChro{
    memory params.mem_req 
    input :
     set val(chro), file(gwas), file(rsinfo) from gwas_format_chro_rs
    output :
      set val(chro), file(gwas),file(outrs) into rsinfo_chro
    script :
      outrs="info_rs_"+chro+".rs"
    """
    zcat $rsinfo | extractrsid_bypos.py --file_chrbp $gwas --out_file $outrs --ref_file stdin --chr $chro --chro_ps ${params.poshead_chro_inforef} --bp_ps ${params.poshead_bp_inforef} --rs_ps ${params.poshead_rs_inforef} --a1_ps ${params.poshead_a1_inforef}  --a2_ps ${params.poshead_a2_inforef} 
    """ 
}

if(params.input_dir!="" || params.input_pat!=''){
print("used plink file")
bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()


}else{
bed=file('bed')
bim=file('bim')
fam=file('fam')
}

fileplk= Channel.create()
Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [a[0], a[1], a[2]] }
    .set { fileplk }

rsinfo_chroplk=rsinfo_chro.combine(fileplk)
process MergeRsGwasChro{
    memory params.mem_req 
    input :
      set val(chro), file(gwas),file(chrors), file(bed), file(bim), file(fam) from rsinfo_chroplk
    output :
      file(outmerge) into gwas_rsmerge
    script :
     outmerge="merge_"+chro+".gwas"
     bfileopt= (params.input_pat!="" || params.input_dir!="") ?  " --bfile "+bed.baseName : ""
     Nheadopt=(params.head_N!="") ? " --N_head ${params.head_N} " : ""
     Freqheadopt=(params.head_freq!="") ? " --freq_head ${params.head_freq} " : ""

     NheadNewopt=(params.headnew_N!="") ? " --Nnew_head ${params.headnew_N} " : ""
     FreqNewheadopt=(params.headnew_freq!="") ? " --freqnew_head ${params.headnew_freq} " : ""
     """
     mergeforrs.py --input_gwas $gwas --input_rs $chrors  --out_file $outmerge --chro_head  ${params.headnew_chr} --bp_head  ${params.headnew_bp} --rs_head ${params.headnew_rs} --chro $chro $bfileopt  $Nheadopt $Freqheadopt $NheadNewopt $FreqNewheadopt  --a1_head ${params.headnew_A1} --a2_head  ${params.headnew_A2}
     """

}
gwas_rsmerge_all=gwas_rsmerge.collect()

process MergeAll{
   memory params.mem_req 
   input :
      file(allfile) from gwas_rsmerge_all
   publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
   output :
      file(fileout) 
   script :
     file1=allfile[0]
     listefiles=allfile.join(" ")
     fileout=params.output
     """
     head -1 $file1 > $fileout
     ls $listefiles  | xargs -n 1 tail -n +2 >> $fileout
     """

}


