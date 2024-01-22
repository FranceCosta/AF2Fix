#!/bin/bash

target_proteins=output/DB_query/target_proteins.tsv

while getopts ":i:u:p:h:P:" flag; do
    case $flag in
        i) target_proteins=${OPTARG};;
        u) SQLuser=${OPTARG};; # user for pfam MySQL
        p) SQLpassword=${OPTARG};; # password for pfam MySQL
        h) SQLhost=${OPTARG};; # host for pfam MySQL
        P) SQLport=${OPTARG};; # port for pfam MySQL
    esac
done

domains=$(tail -n +2 $target_proteins | cut -d ' ' -f1 | uniq)

echo domain,length > domain_lengths.csv

for domain in $domains
do
  mysql -h $SQLhost -u $SQLuser -p$SQLpassword -P $SQLport pfam_35_0 --quick -e "SELECT hmm FROM pfamA_HMM WHERE pfamA_acc='$domain'" | grep -oP -i 'LENG  [0-9]{1,}' | echo $domain,$(cut -d ' ' -f3) >> domain_lengths.csv
done
