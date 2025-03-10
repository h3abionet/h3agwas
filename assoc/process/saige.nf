include { strmem } from "../../modules/fct_groovy.nf"
process COFACTORS_TYPE {
    executor 'local'
    input:
    val args_ch
    val infoargs
    val regenie

    output:
    val CofactStr

    exec :
    def splargs = args_ch.split(",")
    def splinfoargs = infoargs.split(",")
    
    if (infoargs == '') {
        CofactStr = "--covarColList ${args_ch}"
        return CofactStr
    }

    if (splargs.size() != splinfoargs.size()) {
        throw new IllegalArgumentException("Covariates and covariate types are not the same size: ${args_ch} ${infoargs}")
    }

    def Cofactqual = []
    def Cofactquant = []
    def allcov = []

    for (int i = 0; i < splargs.size(); i++) {
        allcov << splargs[i]

        if (splinfoargs[i] == '1') {
            Cofactqual << splargs[i]
        } else if (splinfoargs[i] == '0') {
            Cofactquant << splargs[i]
        } else {
            throw new IllegalArgumentException("Unknown type for ${splargs[i]}: ${splinfoargs[i]}. Use 1 for qualitative, 0 for quantitative.")
        }
    }

    Cofactqual = Cofactqual.join(",")
    Cofactquant = Cofactquant.join(",")
    allcov = allcov.join(",")

    if (regenie == 1) {
        CofactStr = ""
        if (Cofactqual) CofactStr += " --catCovarList ${Cofactqual}"
        if (allcov) CofactStr += " --covarColList ${allcov}"
    } else {
        CofactStr = ""
        if (Cofactqual) CofactStr += " --qCovarColList ${Cofactqual}"
        CofactStr += " --covarColList ${allcov}"
    }
    return CofactStr
}
   def cofactors_type(args_ch,infoargs, regenie) {                                  
      infoargs=""+infoargs                                                      
      splargs=args.split(",")                                                   
      splinfoargs=infoargs.split(",")                                           
      if(infoargs==''){                                                         
        cov="--covarColList "+args                                              
        return(cov)                                                             
      }                                                                         
      if(splargs.size() != splinfoargs.size()){                                 
         System.err.println("covariates and covariate type are not same size : "+args+" "+infoargs)
         System.exit(-11)                                                       
      }                                                                         
      Cofactqual=""                                                             
      Cofactquant=""                                                            
      allcov=""                                                                 
      for (i = 0; i <splargs.size(); i++) {                                     
          /*0 : for quantitatif */                                              
          /* 1 for qualitatif*/                                                 
          if(allcov=="")allcov=splargs[i]                                       
          else allcov+=","+splargs[i]                                           
          if     (splinfoargs[i]=='1'){                                         
              if(Cofactqual=="")Cofactqual=splargs[i]                           
              else Cofactqual+=","+splargs[i]                                   
                                                                                
          }else {                                                               
           if(splinfoargs[i]=='0'){                                             
              if(Cofactquant=="")Cofactquant=splargs[i]                         
              else Cofactquant+=","+splargs[i]                                  
           }else{                                                               
             System.err.println("type args for "+splargs[i]+" doesn't know "+ splinfoargs[i]+"\n 1 for qualitatif arguments\n 0 for quantitatif arguments")
             System.exit(-10)                                                   
          }                                                                     
        }                                                                       
      }                                                                         
     CofactStr=""                                                               
     if(regenie==1){                                                            
       if(Cofactqual!="")CofactStr=" --catCovarList "+Cofactqual                
        if(Cofactquant!="") CofactStr+=" --covarColList "+allcov                
     }else{                                                                     
       if(Cofactqual!="")CofactStr=" --qCovarColList  "+Cofactqual              
        CofactStr+=" --covarColList "+allcov                                    
     }                                                                          
      return(CofactStr)      
}

process bgen_formatsample {                                                 
   label 'R'                                                               
   input :                                                                 
         path(data) 
         path(bgen_sample) 
   output :                                                                 
          path(bgen_sample2)  
          path(bgen_samplesaige2) 
   script :                                                                 
        bgen_sample2=bgen_sample+".modif"                                    
        bgen_samplesaige2=bgen_sample+"_saige.modif"                         
        """                                                                  
        format_samplebgen.r --sample $bgen_sample  --out_sample $bgen_sample2 --out_samplesaige $bgen_samplesaige2 --data $data
        """                                                                  
}     

process checkidd_saige_vcf{
    label 'R'
   memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
   cpus params.max_cpus
    input :
      path(data)
      val(pheno)
      val(pheno_bin)
      path(vcf) 
      tuple path(bed), path(bim), path(fam) 
    output :
      tuple path("${bfileupdate}.bed"), path("${bfileupdate}.bim"), path("${bfileupdate}.fam"), emit : plink
      path(out), emit : pheno
    script :
      out=data.baseName+'_saige.ind'
      bfile=bed.baseName
      bfileupdate=bfile+'_idsaige'
      """
      zcat $vcf |head -1000 |grep "#"| tail -1|awk '{for(Cmt=10;Cmt<=NF;Cmt++)print \$Cmt}' > fileind
      awk '{print \$1"\t"\$1}' fileind > keep
      format_saige_pheno.r --data $data --ind_vcf fileind --out ${out} --pheno ${pheno}  --pheno_bin ${pheno_bin}
      awk '{if( \$1 ~ /^[0-9]+\$/)print \$1"\t"\$4"\t"\$4"\t"\$4}' $bim > keep.range
      plink -bfile  $bfile --keep-allele-order -out $bfileupdate --make-bed --update-ids  ${out}"_updateid" --keep keep --allow-extra-chr --extract range keep.range  --threads ${params.max_cpus}

      """
}

process checkidd_saige{
    label 'R'
   memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
   cpus params.max_cpus
    input :
      path(data)
      val(pheno)
      val(pheno_bin)
      path(bgensample)
      tuple path(bed), path(bim), path(fam)
    output :
      tuple path("${bfileupdate}.bed"), path("${bfileupdate}.bim"), path("${bfileupdate}.fam"), emit :plink
      tuple path(covariates_form),path("${bfileupdate}.fam"), emit : pheno
    script :
      covariates_form=data.baseName+'_saige.ind'
      bfile=bed.baseName
      bfileupdate=bfile+'_idsaige'
      indnbgen=(params.bgen!='') ? "--ind_bgen $bgensample" : ""
      """
      format_saige_pheno.r --data $data $indnbgen --out ${covariates_form}  --pheno ${pheno}  --pheno_bin ${pheno_bin}
      awk '{if( \$1 ~ /^[0-9]+\$/)print \$1"\t"\$4"\t"\$4"\t"\$4}' $bim > keep.range
      plink -bfile  $bfile --keep-allele-order -out $bfileupdate --make-bed --update-ids  ${covariates_form}"_updateid" --extract range keep.range --threads ${params.max_cpus}
      """
}
process  getSaigePheno {
    input:
      tuple path(covariates), path(fam) from data_ch_saige_form
    output:
      file(phef) into saige_data_format_ch
      stdout into pheno_cols_ch_saige
    script:
      covoption = (params.covariates=="") ? "" : " --cov_list ${params.covariates}"
      base = fam.baseName
      phef = "${base}_saige_n.phe"
      """
      all_covariate.py --data  $covariates --inp_fam  $fam $covoption \
                          --pheno ${params.pheno} --phe_out ${phef}  --form_out 2 --nona 1
      """
}


process saige_computed_variance{ 
   memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
   cpus params.max_cpus
   label 'saige'
   input :
      tuple path(data),val(pheno),val(pheno_bin),val(covariable_saige),path(bed), path(bim), path(fam) 
   publishDir "${params.output_dir}/saige/varexp/", overwrite:true, mode:'copy'
   output :
     tuple val(our_pheno) ,path("${output}.rda"), path("${output}.varianceRatio.txt"), val(plk), emit : var
     path("*.report")
   script :
     binpheno = (pheno_bin==1) ? " --traitType=binary " : " --traitType=quantitative"
     Loco = (params.loco==1) ? " --LOCO=TRUE " : " --LOCO=FALSE "
     plk=bed.baseName
     our_pheno    = pheno.replaceAll(/|\/np.\w+/,"").replaceAll(/[0-9]+@@@/,"")
     output=our_pheno+"_var"
     covoption= (covariable_saige=="") ? "" : "  ${covariable_saige} "
     """
     ${params.saige_bin_fitmodel} \
        --plinkFile=$plk \
        --phenoFile=$data \
        --phenoCol=$our_pheno $covoption \
        --sampleIDColinphenoFile=IID \
        --outputPrefix=./$output \
        --nThreads=${params.max_cpus} $Loco \
        --minMAFforGRM=${params.snp_maf_rel} $binpheno \
        ${params.saige_otheropt_step1} ${params.saige_otheropt}  --maxMissing ${params.cut_geno}
      cp .command.sh "${our_pheno}"_saige_step1.cmd.report
      cp .command.log "${our_pheno}"_saige_step1.log.report
      cp .command.err "${our_pheno}"_saige_step1.err.report
     """
  }
/*
process buildindex {
   memory { strmem(params.other_process_mem_req) + 5.GB * (task.attempt -1) }               
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }         
   maxRetries 3                                                           
    label 'utils'
    input :
       file(vcf) from listvcf_ch
    output :
        set file(vcf), file(vcfindex) into listvcf_ind_ch
    script :
       vcfindex=vcf+".csi"
       """
       hostname
       tabix -C -p vcf $vcf
       """
}
*/
process doSaige{
    memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
    cpus params.max_cpus
    label 'saige'
    input :
       tuple val(our_pheno), path(rda), path(varRatio), val(base), path(vcf),path(vcfindex)
   publishDir "${params.output_dir}/saige/log/", overwrite:false, mode:'copy', pattern: "*.report"
    output :
      tuple val(our_pheno),path("$output"), val(base)
      path("*.report")
    script :
     output=vcf.baseName+".res"
     Loco = (params.loco==1) ? " --LOCO=TRUE " : " --LOCO=FALSE "
     imputed=(params.saige_imputed_data==1) ? "--is_imputed_data=TRUE " : '--is_imputed_data=FALSE' 
     moredetail=" --is_output_moreDetails=TRUE " 
     """
      Chro=`zcat $vcf|grep -v "#"|head -1|awk '{print \$1}'`
      ${params.saige_bin_spatest} \
        --vcfFile=$vcf\
        --vcfFileIndex=$vcfindex \
        --vcfField=${params.vcf_field} \
        --minMAF=${params.cut_maf}\
        --minMAC=${params.vcf_minmac} \
        --chrom=\$Chro \
        --GMMATmodelFile=$rda \
        --varianceRatioFile=$varRatio \
        --SAIGEOutputFile=$output $Loco $imputed  $moredetail ${params.saige_otheropt_step2}  ${params.saige_otheropt}    --maxMissing ${params.cut_geno}
      cp .command.sh "${our_pheno}"_saige_step2.cmd.report
      cp .command.log "${our_pheno}"_saige_step2.log.report
      cp .command.err "${our_pheno}"_saige_step2.err.report
     """
}
process doSaigeBgen{
     memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
     cpus params.max_cpus
    label 'saige'
    input :
       tuple val(pheno),path(rda), path(varRatio), val(base),val(Chro),path(bgen), path(bgenindex), path(bgensample)
   publishDir "${params.output_dir}/saige/log/", overwrite:false, mode:'copy', pattern: "*.report"
    output :
      tuple val(our_pheno), path("$output"), val(base)
      path("*.report")
    script :
     output=bgen.baseName+"_"+Chro+".saige"
     Chro=Chro.replace(' ','')
     Chro=(Chro=="23") ? "X": "$Chro"
     Loco = (params.loco==1) ? " --LOCO=TRUE " : " --LOCO=FALSE "
     imputed=(params.saige_imputed_data==1) ? "--is_imputed_data=TRUE --impute_method=${params.saige_impute_method} " : '--is_imputed_data=FALSE' 
     moredetail=" --is_output_moreDetails=TRUE " 
     """
      ${params.saige_bin_spatest} \
        --bgenFile=$bgen    \
        --bgenFileIndex=$bgenindex \
        --sampleFile=$bgensample \
        --AlleleOrder=ref-first \
        --chrom=$Chro \
        --minMAF=${params.cut_maf}\
        --minMAC=${params.vcf_minmac} \
        --GMMATmodelFile=$rda \
	--SAIGEOutputFile=$output $imputed $Loco $moredetail  ${params.saige_otheropt_step2}   ${params.saige_otheropt}  --maxMissing ${params.cut_geno}
      cp .command.sh "${our_pheno}"_saige_step2.cmd.report
      cp .command.log "${our_pheno}"_saige_step2.log.report
      cp .command.err "${our_pheno}"_saige_step2.err.report
     """
}

process getchrobgen{
     label 'utils'
     memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
     cpus params.num_cores

      input :
       tuple path(bgen), path(bgenindex) from bgen_ch_saige
      output :
       tuple env(Chro), path(bgen), path(bgenindex) into bgen_ch_saige_chro
       """
       Chro=`bgenix -g $bgen -list |sed '1,2d'| awk '{print \$3}' |head -1`
       """

}
process doSaigeListBgen{
    memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
    cpus params.num_cores
    label 'saige'
    input :
       tuple val(Chro),path(bgen), path(bgenindex), path(bgensample),val(our_pheno),path(rda), path(varRatio), val(base) from bgen_andvariance_ch
   publishDir "${params.output_dir}/saige/log/", overwrite:false, mode:'copy', pattern: "*.report"
    output :
      set val(our_pheno),file("$output"), val(base) into ch_saige_bychro
      path("*.report")
    script :
     output=bgen.baseName+"_"+Chro+".saige"
     //bin_option_saige= (params.pheno_bin==1) ? " --IsOutputAFinCaseCtrl=TRUE  --IsOutputNinCaseCtrl=TRUE --IsOutputHetHomCountsinCaseCtrl=TRUE " : ""
     Loco = (params.saige_loco==1) ? " --LOCO=TRUE " : " --LOCO=FALSE "
     imputed=(params.saige_imputed_data==1) ? "--is_imputed_data=TRUE --impute_method=${params.saige_impute_method}" : '--is_imputed_data=FALSE' 
     moredetail=(params.pheno_bin==1) ? " --is_output_moreDetails=TRUE " : ""
     """
      ${params.saige_bin_spatest} \
        --bgenFile=$bgen    \
        --bgenFileIndex=$bgenindex \
        --sampleFile=$bgensample \
        --chrom=$Chro \
        --AlleleOrder=ref-first \
        --minMAF=${params.cut_maf}\
        --minMAC=${params.vcf_minmac} \
        --GMMATmodelFile=$rda \
        --varianceRatioFile=$varRatio \
        --SAIGEOutputFile=$output ${Loco} $imputed  $moredetail ${params.saige_otheropt_step2} ${params.saige_otheropt}
      cp .command.sh "${our_pheno}"_saige_step2.cmd.report
      cp .command.log "${our_pheno}"_saige_step2.log.report
      cp .command.err "${our_pheno}"_saige_step2.err.report
     """
}


   

process doSaigePlink{
    memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
    cpus params.max_cpus
    label 'saige'
    input :
       tuple file(bed), file(bim), file(fam),val(our_pheno) , path(rda), path(varRatio), val(base) from plink_andvariance_ch
    output :
      tuple val(our_pheno),file("$output"), val(base) into ch_saige_bychro
      path("*.report")
   publishDir "${params.output_dir}/saige/log/", overwrite:false, mode:'copy', pattern: "*.report"
    script :
     output=bed.baseName+'_'+Chro+".saige"
     Loco = (params.saige_loco==1) ? " --LOCO=TRUE " : " --LOCO=FALSE "
     //bin_option_saige= (params.pheno_bin==1) ? " --IsOutputAFinCaseCtrl=TRUE  --IsOutputNinCaseCtrl=TRUE --IsOutputHetHomCountsinCaseCtrl=TRUE " : ""
     imputed=(params.saige_imputed_data==1) ? "--is_imputed_data=TRUE --impute_method=${params.saige_impute_method} " : '--is_imputed_data=FALSE' 
     moredetail=(params.pheno_bin==1) ? " --is_output_moreDetails=TRUE " : ""
     """

      ${params.saige_bin_spatest} \
        --bedFile=$bed \
        --bimFile=$bim \
        --famFile=$fam \
        --chrom=$Chro \
        --AlleleOrder=alt-first \
        --minMAF=${params.cut_maf}\
        --minMAC=${params.vcf_minmac} \
        --GMMATmodelFile=$rda \
        --varianceRatioFile=$varRatio \
        --SAIGEOutputFile=$output $Loco $imputed  $moredetail ${params.saige_otheropt_step2}  ${params.saige_otheropt}
      cp .command.sh "${our_pheno}"_saige_step2.cmd.report
      cp .command.log "${our_pheno}"_saige_step2.log.report
      cp .command.err "${our_pheno}"_saige_step2.err.report
     """
   }


process join2channel {
  input :
   val(a)
   val(b)
  output :
    tuple val(a),val(b)
  script :
     """ 
     echo 'tmp'   
     """
}

workflow saige{
  take :
   data
   pheno
   pheno_bin
   covariates
   covariates_type
   plink
   plink_rel
   vcf 
   bgen
   bgen_sample
   listchro
 main :
 COFACTORS_TYPE(covariates, covariates_type,channel.of(0))
 if(params.vcf!='' || params.vcf_list!=''){
    checkidd_saige_vcf(data, pheno,pheno_bin,vcf.first(),plink_rel)
    plkinp=checkidd_saige_vcf.out.plink
 }else if(params.bgen || params.bgen_list){
   checkidd_saige(data, pheno, bgen.first(), plink_rel)
   plkinp=checkidd_saige.out.plink
 }else{
   checkidd_saige(data, pheno, channel.fromPath("00"), plink_rel)
   plkinp=checkidd_saige.out.plink
 }
 npheno=params.pheno.split(',').size()
 phenol = pheno.flatMap { list -> list.split(',') }  // Groovy-style lambda for splitting
 pheno_binl = pheno_bin.map { list -> list.split(',') }.flatMap { it.size() == 1 ? it.collect { it } * npheno : it }
 combined = join2channel(phenol,pheno_binl)
 saigeinp=checkidd_saige_vcf.out.pheno.combine(combined).combine(COFACTORS_TYPE.out).combine(plkinp)
 saige_computed_variance(saigeinp)

 if(params.vcf!='' || params.vcf_list!=''){                                     
       doSaige(saige_computed_variance.out.var.combine(vcf))
 }else if(params.bgen!='' || params.bgen_list!=''){
        bgen_formatsample(checkidd_saige.pheno,)
        getchrobgen(bgen)
        doSaigeBgen(saige_computed_variance.out.var.combine(getchrobgen.out))
 }else {
        doSaigePlink(saige_computed_variance.out.var.combine(plink))
 }


}
