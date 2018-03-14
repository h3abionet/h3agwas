#!/bin/bash

hostname


plink --bfile $base --threads 1 --linear $perm $adjust --out $base
