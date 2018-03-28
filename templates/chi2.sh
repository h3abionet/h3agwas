#!/bin/bash

hostname

plink --bfile $base --threads 1 $covariate --assoc $perm $adjust --out $base