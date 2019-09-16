<img src="../aux/H3ABioNetlogo2.jpg"/>

# Pipeline : format data

## Script : `vcf_in_plink.nf`

### objectifs
convert vcf by chromosome, merge  in plink format, added cm  

### arguments :
   * `file_listvcf` : file contains each vcf files to merge, one by line [default : none]
   * `min_scoreinfo` : what score info minimum do you want : [default :0.6]
   * `output_pat` : pattern of output for bed final file [default : out]
   * `output_dir` : directory of output : [default : plink] 
   * `genetic_maps` : genetics maps to added map in bim file, if not provided, map doesn't added in bim, must be not compressed :

```
chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)
1 55550 0 0
1 568322 0 0
1 568527 0 0
```
    * file for [hg19](https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg19.txt.gz) 



