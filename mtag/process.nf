#!/usr/bin/env nextflow

process GetSnpList{                                                           
   time params.big_time                                                         
   input :                                                                      
      tuple path(gwas), val(headchr), val(headbp), val(heada1), val(heada2), val(headrs)
   output :                                                                     
     path("$out"), emit: list_rs_format1 
     path("${out}2"), emit : list_rs_format_2
   script :                                                                     
     out=gwas+".list"                                                           
     chrbp=(headchr!="" && headbp!="") ? ",Pos:${headbp},Chro:${headchr}" : ""
     """                                                                        
     ma_extract_rsid.py --input_file $gwas --out_file $out --info_file rsID:${headrs},A1:${heada1},A2:${heada2}${chrbp} --ldsc
     """                                                                        
} 

 process mergelistsnp {                                                         
    input :                                                                     
      path(allfile) 
    output :                                                                    
      path(mergeall)
    script :                                                                    
      allfile_m=allfile.join(" ")                                               
      mergeall="merge_rs"                                                       
      """                                                                       
      head -1 ${allfile[0]} > $mergeall                                         
      cat $allfile_m|grep -v "SNP"|sort|uniq >> $mergeall                       
      """                                                                       
                                                                                
 } 


                                                                                
process extract_rs_formldsc{                                                    
     input :                                                                    
      tuple path(ldscore), path(posinfo)
     output :                                                                   
       path(out)
     script :                                                                   
       out=ldscore+".rs"                                                        
     """                                                                        
     extract_rs_ldscdb.py --ldsc_input $ldscore --gwas_input $posinfo --out $out
     """                                                                        
 }                                                                              
                                                                                
process merge_allchroscore{                                                     
  input :                                                                       
    path(allfile) 
  output :                                                                      
      path(filers) 
  script :                                                                      
  filers="allrs.update"                                                         
   """                                                                          
   cat ${allfile.join(" ")} >>$filers                                           
   """                                                                          
}         
process format_summary{
 memory params.high_memory
 input :
///dataC/alcohol/v2/multitrait/SI_meta.tsv, rs, p, af, a1, a2, a2, se, beta, chr, bp, /dataC/alcohol/v2/multitrait/work/d7/4ef10125630529ea6a944c2954833a/allrs.update, /dataC/alcohol/v2/multitrait/01, /dataC/alcohol/v2/multitrait/02, /dataC/alcohol/v2/multitrait/03]

   tuple path(filegwas), val(head_rs),val(head_pval), val(head_freq), val(head_A1), val(head_A2), val(head_se), val(head_beta),val(head_chr),val(head_bp),val(head_N),path(file_updaters),path(bed), path(bim), path(fam) 
 output :
   path(out)
 script :
    out=filegwas.toString().replace('-','_')+".format"
    headfreq=head_freq!="" ? " --freq_header  ${head_freq} " : ""
    headN=head_N!="" ? " --n_head ${head_N} " : ""
    plk=bed.baseName
    cmdplk=(bed ==~ /bed/ ? " --bfile $plk " : "")
     """
     gcta_format.py --inp_asso $filegwas  --rs_header ${head_rs} --pval_header ${head_pval} $headfreq --a1_header ${head_A1} --a2_header ${head_A2} --se_header ${head_se} --beta_header ${head_beta} --chro_header ${head_chr}  --out $out --threads ${params.max_plink_cores} --bp_header ${head_bp} $headN $cmdplk --print_pos 1 --file_updaters $file_updaters 
     """
}

///*Mtag*/
process perform_mtag{
   label 'mtag'
   memory params.high_memory
   time params.big_time
 input :
   path(listfile)
   tuple val(headn), val(listn), path(dirwld)
 output :
   tuple val(fnames),path("${out}_trait*")
 script :
   fnames = listfile.join(",")
   out=''
   Ninfo=headn!="" ? " --n_name N " : " --n_value listn" 
   ////--info_min $info_min
   """
   ${params.bin_mtag} --sumstats $fnames --out ./$out --snp_name SNP --beta_name b --se_name se --eaf_name freq --maf_min ${params.cut_maf} --a1_name A1 --a2_name A2 --chr_name chro --p_name p  --bpos_name bp --incld_ambig_snps   ${params.opt_mtag} --force
   """
}

process renames_mtag{
   input :
     tuple val(listi), path(listf), val(outputdir)
   publishDir "${outputdir}", mode:'copy'
   output :
      file("mtag_res/*") 
   script :
      listf2=listf.join(',')
      """
      mkdir -p mtag_res/
      rename_mtag.py --listi $listi --listf $listf2 --dir_out mtag_res
      """
}
//
//mtag_format_ch_fm=mtag_format_ch.flatMap()
////SNP     CHR     BP      A1      A2      Z       N       FRQ     mtag_beta       mtag_se mtag_z  mtag_pval
//
process showMtag{
    memory params.high_memory
    publishDir params.output_dir, overwrite:true, mode:'copy'
    input:
      file(filegwas) 
    output:
      file("${out}*") 
    script:
      out = filegwas.baseName.replaceAll(/_/,'-').replaceAll(/\./,'-')
      """
      general_man.py  --inp $filegwas --phenoname $out --out ${out} --chro_header CHR --pos_header BP --rs_header SNP --pval_header mtag_pval --beta_header mtag_beta --info_prog "Mtag : all file"
      """
}
//
//if(listfilegwas.size()>2){
//list_file2by2=[]
//nbfile=listfilegwas.size()
//for (i = 0; i <(nbfile-1); i++){
//   for( j = i+1; j<nbfile; j++){
//      list_file2by2.add([i, j])
//}
//}
//
//process doMTAG2by2{
//   label 'mtag'
//   memory params.mtag_mem_req
//   time params.big_time
//   input :
//     file(listfile) from listgwasform_col2
//     each poss from list_file2by2
//     publishDir "${params.output_dir}/mtag_2by_2", overwrite:true, mode:'copy'
//     output:
//       file("$output"+"*")
//     script :
//        file1=listfile[poss[0]]
//        file2=listfile[poss[1]]
//        output=""+file1+"_"+file2
//        """
//        ${params.bin_mtag} --sumstats $file1,$file2 --out ./$output --snp_name SNP --beta_name b --se_name se --eaf_name freq --maf_min ${params.cut_maf} --a1_name A1 --a2_name A2 --chr_name chro --p_name p --bpos_name bp --incld_ambig_snps ${params.opt_mtag} --force --z_name z
//        """
//}
//}
//
//def getres(x) {
//  def  command1 = "$x"
//  def  command2 = "head -n 1"
//  def proc1 = command1.execute()
//  def proc2 = command2.execute()
//  def proc = proc1 | proc2
//  proc.waitFor()
//  res ="${proc.in.text}"
//  return res.trim()
//}
//
//nextflowversion =getres("nextflow -v")
//if (workflow.repository)
//  wflowversion="${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
//else
//  wflowversion="A local copy of the workflow was used"
//
//report_ch = report_mtag.flatten().toList()
//
//process doReport {
//  label 'latex'
//  input:
//    file(reports) from report_ch
//  publishDir params.output_dir, overwrite:true, mode:'copy'
//  output:
//    file("${out}.pdf")
//  script:
//    out = params.output+"-report"
//    these_phenos     = "mtag"
//    these_covariates = ""
//    config = getConfig()
//    images = workflow.container
//    texf   = "${out}.tex"
//    template "make_assoc_report.py"
//
//}
//
