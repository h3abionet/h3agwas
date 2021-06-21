pi_hat          = params.pi_hat
super_pi_hat    = params.super_pi_hat
cut_diff_miss   = params.cut_diff_miss
f_lo_male       = params.f_lo_male
f_hi_female     = params.f_hi_female

plink_mem_req = params.plink_mem_reqq 
other_mem_req = params.other_mem_reqq
max_plink_cores = params.max_plink_coress

sexinfo = ""
extrasexinfo = "--must-have-sex"

K = "--keep-allele-order"

process getListOfDuplicatePositions {

    label 'plink'

    input:
        tuple path(bed), path(bim), path(fam)

    output:
        tuple path('plink.dupvar'), path('plink.hh')

    script:
        """
        plink \
            --bfile ${params.cohortName} \
            --list-duplicate-vars \
            ids-only \
            suppress-first
        """
}

process removeSamplesWithPoorClinicalData {
    label 'plink'
    
    input:
        tuple path(bed), path(bim), path(fam), path(poorSamples)
    
    output:
        path "${output}.{bed,bim,fam}"        
    
    script:
        output = "${params.cohortName}.sampleFiltered"
        """
        plink \
            --bfile ${params.cohortName} \
            --remove ${poorSamples} \
            --make-bed \
            --out ${output}
        """
}

process removeDuplicatedSnvPositions {

    label 'plink'
    publishDir "${params.outputDir}/quality-control", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)
    tuple path('plink.dupvar'), path('plink.hh')

    output:
        path("${output}.{bed,bim,fam,log}")

    script:
        output = "${params.cohortName}"
        """
        plink \
            --bfile ${params.cohortName}.sampleFiltered \
            --make-bed \
            -exclude plink.dupvar \
            --out ${output}
        """
}

process identifyIndivDiscSexinfo {

  memory other_mem_req
  publishDir "${params.output_dir}/IndivDiscSexinfo", overwrite:true, mode:'copy'

  input:
     tuple path(bed), path(bim), path(fam),path(log)

  output:
     tuple path(lmiss), path(imiss), path(sexcheck_report), path("${base}.hwe"), path(logfile)

  script:
    base = bed.baseName
    logfile= "${base}.badsex"
    sexcheck_report = "${base}.sexcheck"
    imiss  = "${base}.imiss"
    lmiss  = "${base}.lmiss"
    """
       plink --keep-allele-order --bfile $base --hardy --check-sex $f_hi_female $f_lo_male --missing  --out $base
       head -n 1 ${base}.sexcheck > $logfile
       grep  'PROBLEM' ${base}.sexcheck >> $logfile
    """
}

process generateSnpMissingnessPlot {
  tag "Missingness plot"

  memory other_mem_req
  publishDir "${params.output_dir}/miss", overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
     tuple path(lmissf), path(imiss), path(sexcheck_report), path(hwe), path(logfile)

  output:
     path output

  script:
    input  = lmissf
    base   = lmissf.baseName
    label  = "SNPs"
    output = "${base}-snpmiss_plot".replace(".","_")+".pdf"
    template "missPlot.py"
}

process generateIndivMissingnessPlot {

  memory other_mem_req
  publishDir "${params.output_dir}/miss", overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
     tuple path(lmissf), path(imissf), path(sexcheck_report), path(hwe), path(logfile)
  
  output:
    path output
  
  script:
    input  = imissf
    base   = imissf.baseName
    label  = "samples"
    output = "${base}-indmiss_plot".replace(".","_")+".pdf"
    template "missPlot.py"
}

process getInitMAF {

  memory plink_mem_req

  input:
    tuple path(bed), path(bim), path(fam)

  output:
     path "${newbase}.frq"

  script:
    base = bed.baseName
    newbase = base.replace(".","_")

    """
    plink --bfile $base --freq --out $newbase
    """
}

process showInitMAF {

  memory other_mem_req
  publishDir "${params.output_dir}/InitMAF", overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
     path freq

  output:
     tuple path("${base}.pdf"), path("${base}.tex")

  script:
    base = freq.baseName+"-initmaf"
    base = base.replace(".","_")
    template "showmaf.py"
}

process showHWEStats {

  memory other_mem_req
  publishDir "${params.output_dir}/HWEStats", overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
     tuple path(lmissf), path(imissf), path(sexcheck_report), path(hwe), path(logfile)

  output:
     tuple path("${base}.pdf"), path("${base}-qq.pdf"), path("${base}.tex")

  script:
    base = hwe.baseName+"-inithwe"
    base = base.replace(".","_")
    template "showhwe.py"
}

process removeQCPhase1 {
  tag "Removing QC Phase 1"

  memory plink_mem_req
  publishDir "${params.output_dir}/QCPhase1", overwrite:true, mode:'copy'

  input:
    tuple path(bed), path(bim), path(fam)

  output:
    //tuple path("${output}*.bed"), path("${output}*.bim"), path("${output}*.fam"), path("qc1.out"), path("${output}.irem")
    tuple path("${output}*.bed"), path("${output}*.bim"), path("${output}*.fam"), path("${output}.irem")
  
  script:
     base = bed.baseName
     output = "${base}-c".replace(".","_")
     """
     # remove really realy bad SNPs and really bad individuals
     plink $K --autosome --bfile $base $sexinfo --mind 0.1 --geno 0.1 --make-bed --out temp1
     plink $K --bfile temp1  $sexinfo --mind $params.cut_mind --make-bed --out temp2
     plink $K --bfile temp2  $sexinfo --geno $params.cut_geno --make-bed --out temp3
     plink $K --bfile temp3  $sexinfo --maf $params.cut_maf --make-bed --out temp4
     plink $K --bfile temp4  $sexinfo --hwe $params.cut_hwe --make-bed  --out $output 
     cat *log > logfile
     touch tmp.irem
     cat *.irem > ${output}.irem
     # qc1logextract.py logfile ${output}.irem > qc1.out     
  """
}

process pruneForIBD {

  cpus max_plink_cores
  memory plink_mem_req
  publishDir "${params.output_dir}/pruneForIBD", overwrite:true, mode:'copy'

  input:
    //tuple path(bed), path(bim), path(fam), path(qc1), path(irem)
    tuple path(bed), path(bim), path(fam), path(irem)

  output:
    path "${outf}.genome"

  script:
    base   = bed.baseName
    outf   =  base.replace(".","_")

    if (params.high_ld_regions_fname != "") {
      ldreg_ch=Channel.fromPath(params.high_ld_regions_fname)
      range = " --exclude range $ldreg_ch"
    }
    else{
      range =""
    }

    """
     plink --bfile $base --threads $max_plink_cores --autosome $sexinfo $range --indep-pairwise 60 5 0.2 --out ibd
     plink --bfile $base --threads $max_plink_cores --autosome $sexinfo --extract ibd.prune.in --genome --out ibd_prune
     plink --bfile $base --threads $max_plink_cores --autosome $sexinfo --extract ibd.prune.in --genome --min $pi_hat --out $outf
      echo DONE
     """
}

// run script to find a tuple of individuals we can remove to ensure no relatedness
//  Future - perhaps replaced with Primus
process findRelatedIndiv {

  memory other_mem_req
  publishDir "${params.output_dir}/findRelatedIndiv", overwrite:true, mode:'copy'

  input:
    tuple path(lmiss), path(missing), path(sexcheck_report), path(hwe), path(logfile)
    path(ibd_genome)

  output:
     path outfname

  script:
     base = missing.baseName
     outfname = "${base}-fail_IBD".replace(".","_")+".txt"
     template "removeRelInds.py"
}

process calculateSampleHeterozygosity {

   memory plink_mem_req
   publishDir "${params.output_dir}/Heterozygosity", overwrite:true, mode:'copy'

   input:
      // tuple path(bed), path(bim), path(fam), path(qc1), path(irem)
      tuple path(bed), path(bim), path(fam), path(irem)


   output:
      tuple path("${hetf}.het"), path("${hetf}.imiss")

   script:
      base = bed.baseName
      hetf = "${base}".replace(".","_")
      """
        plink --bfile $base  $sexinfo --het --missing  --out $hetf
      """
}

process generateMissHetPlot {

  memory other_mem_req
  publishDir "${params.output_dir}/MissHetPlot", overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
    tuple path(het), path(imiss)
  
  output:
    path output

  script:
    base = imiss.baseName
    output  = "${base}-imiss-vs-het".replace(".","_")+".pdf"
    template "missHetPlot.py"
}

// Find those who have bad heterozygosity
process getBadIndivsMissingHet {

  memory other_mem_req
  publishDir "${params.output_dir}/BadIndivsMissingHet", overwrite:true, mode:'copy', pattern: "*.txt"

  input:
    tuple path(het), path(imiss)

  output:
    path outfname

  script:
    base = het.baseName
    outfname = "${base}-fail_het".replace(".","_")+".txt"
    template "select_miss_het_qcplink.py"
}

process noSampleSheet {
    tag "Importing samplesheets"

    output:
     path("poorgc10.lst")

    script:
      template 'sampleqc-minimal.py'
}


// process noSampleSheet {
//     tag "Importing samplesheets"

//     output:
//      tuple path("poorgc10.lst"), path("plates")

//     script:
//       """
//       mkdir -p plates
//       sampleqc.py 0 0 0 poorgc10.lst plates/crgc10.tex
//       """
// }

process removeQCIndivs {

  memory plink_mem_req

  input:
    path(f_miss_het)
    path(rel_indivs)
    tuple path(lmissf), path(imiss), path(sexcheck_report), path(hwe), path(f_sex_check_f)
    path(poorgc)
    // tuple path(bed), path(bim), path(fam), path(qc1), path(irem)
    tuple path(bed), path(bim), path(fam), path(irem)
    
  output:
    path "${out}.{bed,bim,fam}"

  script:
   base = bed.baseName
   out  = "${base}-c".replace(".","_")
    """
     cat $f_sex_check_f $rel_indivs $f_miss_het $poorgc | sort -k1 | uniq > failed_inds
     plink $K --bfile $base $sexinfo --remove failed_inds --make-bed --out $out
     mv failed_inds ${out}.irem
  """
}

mperm_header=" CHR                               SNP         EMP1         EMP2 "

col    = params.case_control_col

process calculateSnpSkewStatus {

  memory plink_mem_req
  cpus max_plink_cores

  input:
    tuple path(bed), path(bim), path(fam), path(ch)

  output:
    tuple path("${base}.missing"), path(mperm), path("${base}.hwe")

  script:
   base  = bed.baseName
   out   = base.replace(".","_")
   mperm = "${base}.missing.mperm"
   diffpheno = "--pheno cc.phe --pheno-name $col"
   phe   = ch
   """
    cp $phe cc.phe
    plink --threads ${max_plink_cores} --autosome --bfile $base $sexinfo $diffpheno --test-missing mperm=10000 --hardy --out $out
    if ! [ -e $mperm ]; then
       echo "$mperm_header" > $mperm
    fi
   """
}

process generateDifferentialMissingnessPlot {

   memory other_mem_req
   publishDir "${params.output_dir}/DMP", overwrite:true, mode:'copy', pattern: "*.pdf"

   input:
     tuple path(clean_missing), path(mperm), path(hwe)
   
   output:
      path output

   script:
       input = clean_missing
       base  = clean_missing.baseName.replace(".","_").replace("-nd","")
       output= "${base}-diff-snpmiss_plot.pdf"
       template "diffMiss.py"
 }

 process findSnpExtremeDifferentialMissingness {

  echo true
  memory other_mem_req

  input:
    tuple path(clean_missing_1), path(clean_missing), path(hwe)

  output:
     path failed

  script:

    cut_diff_miss = params.cut_diff_miss
    missing = clean_missing
    base     = missing.baseName.replace("-.*","").replace(".","_")
    probcol = "EMP2"
    failed   = "sampleA-failed_diffmiss.snps"

    """echo $cut_diff_miss and $missing and $base and $probcol and $failed and 1000 """
    template "select_diffmiss_qcplink.py"
}

process removeSkewSnps {

  memory plink_mem_req
  publishDir "${params.output_dir}/skewSnps", overwrite:true, mode:'copy'  

  input:
    tuple path(bed), path(bim), path(fam)
    path(failed)

  output:
    tuple path("${output}.bed"), path("${output}.bim" ), path("${output}.fam"), path("${output}.log")
  
  script:
  base = bed.baseName
  output = params.output.replace(".","_")
  """
  plink $K --bfile $base $sexinfo --exclude $failed --make-bed --out $output
  """
}

process calculateMaf {

  memory plink_mem_req
  publishDir "${params.output_dir}/Maf", overwrite:true, mode:'copy', pattern: "*.frq"

  input:
    tuple path(bed), path(bim), path(fam), path(log)

  output:
    path "${base}.frq"

  script:
    base = bed.baseName
    out  = base.replace(".","_")
    """
      plink --bfile $base $sexinfo  --freq --out $out
    """
}

process generateMafPlot {

  memory other_mem_req
  publishDir "${params.output_dir}/Maf", overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
    path(input)

  output:
    path output

  script:
    base    = input.baseName
    output  = "${base}-maf_plot.pdf"
    template "mafplot.py"
}

process findHWEofSNPs {

  memory other_mem_req

  input:
    tuple path(missing), path(mperm), path(hwe)

  output:
     path output

  script:
    base   = hwe.baseName.replace(".","_")
    output = "${base}-unaff.hwe"
    """
      head -1 $hwe > $output
      grep 'UNAFF' $hwe >> $output
    """
}

process generateHwePlot {

  memory other_mem_req
  publishDir "${params.output_dir}/HwePlot", overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
    path(unaff)
  
  output:
    path(output)

  script:
    input  = unaff
    base   = unaff.baseName.replace(".","_")
    output = "${base}-hwe_plot.pdf"
    template "hweplot.py"
}