
process merge_sumstat{
    input :
       tuple val(out),val(pheno),path(list_file)
    output :
      tuple val(out),val(pheno), path(out) 
    publishDir "${params.output_dir}/assoc/saige/", overwrite:true, mode:'copy' 
    script :
      fnames = list_file.join(" ")
      file1  = list_file[0]
      outf=out+'_'+pheno+'.assoc'
      """
      head -1 $file1 > $out
      tail -n +2 -q $fnames >> $out
      """
 }
