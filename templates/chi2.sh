#!/bin/bash

echo $num_assoc_cores

plink --bfile $base --threads $num_assoc_cores --assoc $fisher $perm $adjust --out $base

