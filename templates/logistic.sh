#!/bin/bash


hostname 
echo $num_assoc_cores

plink --bfile $base $covariate $pheno --threads ${num_assoc_cores} --logistic $perm $adjust --all-pheno --out $base
