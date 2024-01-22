#!/bin/bash

# get all domains found in each protein

target_proteins=output/DB_query/target_proteins.tsv
output=data/protein_info.csv

while getopts ":i:u:p:h:o:P:" flag; do
    case $flag in
        i) target_proteins=${OPTARG};;
        u) SQLuser=${OPTARG};; # user for pfam MySQL database
        p) SQLpassword=${OPTARG};; # password for pfam MySQL database
        h) SQLhost=${OPTARG};; # host for pfam MySQL
        o) output=${OPTARG};;
        P) SQLport=${OPTARG};;
    esac
done

proteins=$(tail -n +2 $target_proteins | cut -d ' ' -f2 | uniq)

echo pfamA_acc,pfamA_id,comment > $output

for protein in $proteins
do
    echo $protein
    mysql -h $SQLhost -u $SQLuser -p$SQLpassword -P $SQLport pfam_35_0 \
    --quick -e "SELECT pfamA_acc,uniprot_acc,seq_start,seq_end FROM uniprot_reg_full WHERE uniprot_acc='$protein'" \
    | tail -n +2 | tr '\t' ',' >> $output
done
