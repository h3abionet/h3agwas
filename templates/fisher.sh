#!/bin/bash

hostname


plink --bfile $base $covariate $pheno_cmd --threads $num_assoc_cores --fisher $perm $adjust --out $base
