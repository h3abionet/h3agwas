include { strmem } from "./fct_groovy.nf"

process bgen_formatsample {
   label 'R'
   input :
         path(data)
         path(bgen_sample)
   output :
          path(bgen_sample2), emit : regenie
          path(bgen_samplesaige2), emit : saige_gcta
   script :
        bgen_sample2=bgen_sample+".modif"
        bgen_samplesaige2=bgen_sample+"_saige.modif"
        """
        format_samplebgen.r --sample $bgen_sample  --out_sample $bgen_sample2 --out_samplesaige $bgen_samplesaige2 --data $data
        """
}

process indexbgen {
   label 'utils'
   input :
     path(bgen) 
   output :
        tuple path(bgen), path("${bgen}.bgi") 
   """
      bgenix -g $bgen -index
   """
}
process getchrobgen{
     label 'utils'
     memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
     input :
       tuple path(bgen), path(bgenindex) 
      output :
       tuple env(Chro), path(bgen), path(bgenindex)
       """
       Chro=`bgenix -g $bgen -list |sed '1,2d'| awk '{print \$3}' |head -1`
       """
}
