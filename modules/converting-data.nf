process illumina2lgen {
    maxForks params.max_forks
    memory params.indiv_memory_req
    time params.time_req

    input:
       tuple file(report), file(array)
    output:
       tuple file("${output}.ped"), file("${output}.map")
    script:
        samplesize = params.samplesize
        idpat = params.idpat
        output = report.baseName
        """
        hostname
        topbottom.py $array $report $samplesize '$idpat' $output
        """
}