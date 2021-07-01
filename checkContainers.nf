/*
 * CHECK CONTAINERS
 * ================
 *
 ******************************************************/
nextflow.enable.dsl=2

workflow {
    checkPerlContainer() | view()
    checkPlinkContainer() | view()
    checkTidyverseContainer() | view()
    checkQqmanContainer() | view()
    checkPlink2Container() | view()
    //checkBcftoolsContainer() | view()
    checkDatatableContainer() | view()
}

process checkPerlContainer {
    label 'perl'
    executor 'local'

    output:
        stdout

    script:
        """
        perl --version &> /dev/null && /usr/bin/printf "perl \\xE2\\x9C\\x94\\n"
        """
}

process checkPlinkContainer {
    label 'plink'
    executor 'local'

    output:
        stdout

    script:
        """
        plink --help &> /dev/null && /usr/bin/printf "plink \\xE2\\x9C\\x94\\n"
        """
}

process checkTidyverseContainer {
    label 'tidyverse'
    executor 'local'

    output:
        stdout

    script:
        """
        R --version &> /dev/null && /usr/bin/printf "tidyverse \\xE2\\x9C\\x94\\n"
        """
}

process checkQqmanContainer {
    label 'qqman'
    executor 'local'

    output:
        stdout

    script:
        """
        R --version &> /dev/null && /usr/bin/printf "qqman \\xE2\\x9C\\x94\\n"
        """
}

process checkPlink2Container {
    label 'plink2'
    executor 'local'

    output:
        stdout

    script:
        """
        plink2 --help &> /dev/null && /usr/bin/printf "plink2 \\xE2\\x9C\\x94\\n"
        """
}

process checkBcftoolsContainer {
    label 'bcftools'
    executor 'local'

    output:
        stdout

    script:
        """
        bcftools --help &> /dev/null && /usr/bin/printf "bcftools \\xE2\\x9C\\x94\\n"
        """
}

process checkDatatableContainer {
    label 'datatable'
    executor 'local'

    output:
        stdout

    script:
        """
        python --version &> /dev/null && /usr/bin/printf "python \\xE2\\x9C\\x94\\n"
        """
}






