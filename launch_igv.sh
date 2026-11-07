#!/bin/bash

data_dir=/ccb/salz4-4/dbass13/expr_data/

ann=${data_dir}gene_specificity/gtex_v11_tx_protein_coding/extreme_tau_deviation.gtf
ref=/ccb/salz4-4/dbass13/recursive_introns_data/GRCh38/hg38.fa

/ccb/sw/bin/igv -g $ref $ann
