#!/bin/bash

echo $num_assoc_cores

plink --bfile $base --threads 1 --logistic $perm $adjust --out $base

