#!/bin/bash

echo $num_assoc_cores

plink --bfile $base --threads $num_assoc_cores --assoc $perm $adjust --out $base

