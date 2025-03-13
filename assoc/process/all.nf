
process check_pheno_bin{
  label 'R'
  input  :
     val(listpheno)
     val(listpheno_bin)
     path(data)
    val(software)
  output :
     path(out), emit : data
     env bin_pheno, emit: bin_pheno
  script :
     out=data+'.checkplk'
     """
     check_pheno_bin.r --data $data --listpheno $listpheno --listpheno_bin ${listpheno_bin} --out $out --software $software
     bin_pheno=`cat bin_pheno`
     """
}
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

process merge_sumstat{
    input :
       tuple val(out),val(pheno),path(list_file)
    output :
      tuple val(out),val(pheno), path(outf) 
    publishDir "${params.output_dir}/assoc/saige/", overwrite:true, mode:'copy' 
    script :
      fnames = list_file.join(" ")
      file1  = list_file[0]
      outf=out+'_'+pheno+'.assoc'
      """
      head -1 $file1 > $outf
      tail -n +2 -q $fnames >> $outf
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
