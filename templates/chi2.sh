#!/bin/bash

hostname

plink --bfile $base --threads ${num_assoc_cores} $covariate $pheno --assoc $perm $adjust --out $base
