include {perform_mtag;GetSnpList;mergelistsnp;merge_allchroscore;extract_rs_formldsc;format_summary;renames_mtag;show_mtag;doReport} from './process.nf'

workflow mtag{
  take :
    outputdir
  main :
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
    perform_mtag(format_summary.out.collect(),channel.of(params.sumstat_head_n).combine(channel.of(params.sumstat_list_n)).combine(Channel.fromPath("${params.dir_ref_ld_chr}/",type:"dir")), channel.of("$outputdir"))
    renames_mtag(perform_mtag.out.sumstat.combine(channel.of("$outputdir/mtag/")))
    show_mtag(renames_mtag.out.flatMap{it -> it}.combine(channel.of("$outputdir/mtag/")))
    doReport(show_mtag.out.flatten().toList(), channel.of("$outputdir/mtag/"), channel.of("${params.output}"))
    sumstat_mtag=channel.from("mtag").combine(format_summary.out)
    emit :
      sumstat_mtag = sumstat_mtag
      sumstat_clean = format_summary.out
}
