#!/bin/bash


plink --bfile $base --assoc fisher $adjust --out $base

