
workflow getparams_rgc{
   take :
     plink 
     gc
   main :
    if(gc==null){
       dl_gwascat(channel.of(params.gwas_cat_ftp), channel.of("gwascat"), channel.of("${params.output_dir}/utils/data/"))
       gc=dl_gwascat.out
    } 
   if(plink==null){

   }
}

workflow replication_gc {


}
