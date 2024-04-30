







include { qc } from "./qc/qc.nf"
include { assoc} from "./assoc/assoc.nf"

workflow {


  main:
     if (params.action == "qc") {
        qc()
      }

}
