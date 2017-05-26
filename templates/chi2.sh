#!/bin/bash

echo $num_assoc_cores

plink --bfile $base --threads 1 --assoc $perm $adjust --out $base

