
println "In QC top"


process proc_test {
    input:
    val(num)

    output:
      stdout
    script:
      """
        double.py $num
      """

}


workflow qc {
    df = Channel.of(1,2,3)

    main:
    proc_test(df)
    


}
