include {merge_allchroscore;GetSnpList;mergelistsnp;extract_rs_formldsc;format_summary} from '../mtag/process.nf'
include {munge_sumstat;doLDSC;DoCorrLDSC} from './process.nf'


workflow ldsc {
 take :
   // format smmary 
  summary
  outputdir
 main :
  if(summary==null){
      listgwas=channel.from(params.file_gwas.split(',')).flatMap{it->file(it)}
      listgwas_mtag_snplist=listgwas.combine(channel.of(params.sumstat_head_chr)).combine(channel.of(params.sumstat_head_bp)).combine(channel.of(params.sumstat_head_a1)).combine(channel.of(params.sumstat_head_a2)).combine(channel.of(params.sumstat_head_rs))
      GetSnpList(listgwas_mtag_snplist)
      mergelistsnp(GetSnpList.out.list_rs_format_2.collect())
      check_rs_ldsc=Channel.fromPath("${params.dir_ref_ld_chr}/*.ldscore.gz", checkIfExists:true).combine(mergelistsnp.out)
      extract_rs_formldsc(check_rs_ldsc)
      merge_allchroscore(extract_rs_formldsc.out.collect())
      //val{head_rs},val(head_pval), val(head_freq), val(head_A1), val(head_A2), val(head_se), val(head_beta),val(head_chr),val(head_bp)
      plink=channel.fromPath("01").combine(channel.fromPath("02")).combine(channel.fromPath("03"))
      gwasinfo=listgwas.combine(channel.of(params.sumstat_head_rs)).combine(channel.of(params.sumstat_head_pval)).combine(channel.of(params.sumstat_head_freq)).combine(channel.of(params.sumstat_head_a1)).combine(channel.of(params.sumstat_head_a2)).combine(channel.of(params.sumstat_head_se)).combine(channel.of(params.sumstat_head_beta)).combine(channel.of(params.sumstat_head_chr)).combine(channel.of(params.sumstat_head_bp)).combine(channel.of(params.sumstat_head_n)).combine(merge_allchroscore.out).combine(plink)
      format_summary(gwasinfo)
      summary=format_summary.out
   }
//   summary.view()
  listn=summary.combine(channel.of(""))
  if(params.sumstat_list_n!='')listn=summary.merge(channel.from(params.sumstat_list_n.split(',')))
  listn.view() 
  if(params.munge_keep=="")mungekeep=channel.fromPath("0")  else mungekeep=channel.fromPath(params.munge_keep,checkIfExists:true)
  listnmunge=listn.combine(mungekeep)
  munge_sumstat(listnmunge.combine(channel.of("$outputdir/sumstat_mumge/")))
  doLDSC(munge_sumstat.out.sumstat.combine(Channel.fromPath("${params.dir_ref_ld_chr}/")).combine(channel.of("$outputdir/ldsc_h2/")))
  DoCorrLDSC(munge_sumstat.out.sumstat.collect(), Channel.fromPath("${params.dir_ref_ld_chr}/"), channel.of("$outputdir/ldsc_h2/"),channel.of("${params.output}_h2r2"))
}
