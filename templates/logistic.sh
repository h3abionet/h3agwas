#!/bin/bash


hostname 
echo $num_assoc_cores

plink --bfile $base $covariate --threads 1 --logistic $perm $adjust --out $base
