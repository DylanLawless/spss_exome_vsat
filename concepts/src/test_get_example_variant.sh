#!/bin/bash

cd ~/web/ACMGuru/data

# get the header line for example data
cat ./study_v1_chr21_42411318_43411317.csv | head -n1 > example_variant.csv

# get the first stop variant in example data
grep "stop" ./study_v1_chr21_42411318_43411317.csv | head -n1 >> example_variant.csv

