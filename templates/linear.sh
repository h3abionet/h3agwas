#!/bin/bash

hostname


plink --bfile $base $covariate $pheno  --linear $perm $adjust --out $base
