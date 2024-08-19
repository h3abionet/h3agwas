
process update_rsgwas{
  input :
    tuple path(gwas), val(N),path(rsupdate)
  output :
    tuple path(newgwas), val(N),path(newfilepos)
  script :
    newgwas=gwas+'.updaters'
    newfilepos=gwas+'.listpos'
    nval=(params.head_n=='') ? "" : " --n_header ${params.head_n}"
    info=(params.head_info=='') ? "" : " --info_header ${params.head_info} "
    """
    updaters_gwasldsc.py --gwas $gwas --rstoupdate $rsupdate --a1_header ${params.head_A1} --a2_header ${params.head_A2} --chro_header ${params.head_chr} --bp_header ${params.head_bp} --out $newgwas  --rs_header ${params.head_rs} --out_pos $newfilepos --beta_header ${params.head_beta} --se_header ${params.head_se} --af_header ${params.head_freq} $nval --p_header ${params.head_pval}  $info
    """
}

 process munge_sumstat{
   label 'py2ldsc'
   memory params.high_memory
   input :
      tuple path(gwas), val(Nind),path(file_mergeallele), val(outputdir)
   publishDir "${outputdir}/",mode:'copy'
   output :
    path("$out"+".sumstats.gz"), emit:sumstat
    path("$out"+".log"), emit: log
   script :
     NInfo=Nind!="" ? " --N $Nind " : "--N-col N "
     out=gwas.baseName.replace('_','-')+"_mg"
     gwasf=gwas.baseName
     mergeallele=""
     if(params.munge_keep!='')mergeallele="--merge-alleles $file_mergeallele"
     """
     ${params.munge_sumstats_bin} --sumstats $gwas $NInfo --out $out --snp SNP --p p --frq freq --info-min ${params.cut_info} --maf-min ${params.cut_maf} --a1 A1 --a2 A2 --ignore z $mergeallele  --chunksize 500000  #--a1-inc true
     """
}

process doLDSC{
   label 'py2ldsc'
   input :
      tuple path(gwasf), path(dir_ld), val(outputdir)
   publishDir "${outputdir}/", mode:'copy'
   output :
     tuple val(out), path("${out}.log")
   script :
     out=gwasf.baseName+"_ldsc"
     """
     ${params.bin_ldsc} --h2 $gwasf --ref-ld-chr $dir_ld/ --w-ld-chr $dir_ld/ --out $out ${params.ldsc_h2opt}
     """
}

 process DoCorrLDSC{
  label 'py2ldsc'
  time params.big_time
  memory params.high_memory
  input :
    path(listfilegwas)
    path(dirredld)
    val(outputdir)
    val(outpat)
    each pos 
  publishDir "${outputdir}/", mode:'copy'
  output :
    file("$outpat"+".log")
  script : 
    file1=listfilegwas[pos[0]]                                                 
    file2=listfilegwas[pos[1]]      
    listfile=file1+','+file2
    outpat = "${file2}_${file1}.cor"
    """
    ${params.bin_ldsc} --rg $listfile --ref-ld-chr $dirredld/ --w-ld-chr  $dirredld/ --out  $outpat ${params.ldsc_h2opt}
    """
  }
