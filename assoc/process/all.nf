
process merge_sumstat{
    input :
       tuple val(out),file(list_file)
    output :
      tuple val(out), path(assoc) 
    script :
      fnames = list_file.join(" ")
      file1  = list_file[0]
      """
      head -1 $file1 > $assoc
      tail -n +2 -q $fnames >> $out
      """
 }
