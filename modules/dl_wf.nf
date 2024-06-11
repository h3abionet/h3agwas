include {dl_fasta} from './dl.nf'
include {checkfasta} from './utils.nf'
include {fileexist_param} from './fct_groovy.nf'


workflow df_fasta_wf {
  take :
    fasta_file
    build
    link
  main :
   if(build=='' && link=='' && fasta_file==''){
    System.exit(2) 
   }
  if(fasta_file!=''){
    fileexist_param(fasta_file,'fasta')
    fasta_tmp=channel.fromPath(fasta_file)
     checkfasta(fasta_tmp,channel.of("${params.output_dir}/utils/data/"))
     fasta=checkfasta.out.fasta
     fasta_index = checkfasta.out.fasta_index
  }else {
   link = params.ftp_fasta
   if(link==''){
    if(build=='hg19'){
      link=params.ftp_fasta_hg19
     }
    if(build=='37'){

      link = params.ftp_fasta_37
    }
    if(build=='38'){
        link=params.ftp_fasta_38
     }
   }
   if(link==''){
    System.exit(2)
   }
  if(build==''){
    build=File(link).name
    out=build
   }else{
      out=build+'.fa.gz'
   }
   dl_fasta(channel.of(build), channel.of(link), channel.of(out),channel.of("${params.output_dir}/utils/data/"))
   fasta=dl_fasta.out.fasta
   fasta_index = dl_fasta.out.fasta_index
 }
 emit :
  fasta = fasta
  fasta_index = fasta_index
}

