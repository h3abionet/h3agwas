
process MD5_plk {
  input:
     path(plink)
  output:
     path(out)
  script:
       bed = plink[0]
       bim = plink[1]
       fam = plink[2]
       out  = "${plink[0].baseName}.md5"
       template "md5_plk.py"
}

