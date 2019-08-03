#!/bin/bash

hostname


plink --bfile $base $covariate $pheno_cmd --threads $num_assoc_cores --linear $perm $adjust --out $base
