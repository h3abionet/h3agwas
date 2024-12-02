
workflow read_sumstat {
/*read sumstat*/
if(params.file_gwas!=''){
    listgwas=channel.from(params.file_gwas.split(',')).flatMap{it->file(it)}
    listgwas_mtag_snplist=listgwas.combine(channel.of(params.sumstat_head_chr)).combine(channel.of(params.sumstat_head_bp)).combine(channel.of(params.sumstat_head_a1)).combine(channel.of(params.sumstat_head_a2)).combine(channel.of(params.sumstat_head_rs)).combine(channel.of(params.sumstat_head_beta)).combine(channel.of(params.sumstat_head_se)).combine(channel.of(params.sumstat_head_p)).combine(channel.of(params.sumstat_head_z))
}


}
