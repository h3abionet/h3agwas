#!/bin/bash

hostname
echo $num_assoc_cores

plink --bfile $base --threads 1 --linear $perm $adjust --out $base

