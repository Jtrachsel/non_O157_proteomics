#!/bin/bash
set -e

mkdir reference/psort_res
psortb -n -r reference/psort_res -o terse -v -i reference/detected_proteins.faa
cp reference/psort_res/20220817180014_psortb_gramneg.txt output/psort_res.tsv

