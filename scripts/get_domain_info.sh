#!/bin/bash

password=""
port=""
user=""
host=""

# Get relevant info associated to every domain 
# conda activate pymol

# get domains with bimodal distribution

#echo "Getting domains with bimodal distribution"
#python3 ../scripts/get_all_data_from_DB.py --number_to_consider 50 --proteins_per_class 10
target_proteins=data/afdb_protein_info.csv
domains=$(tail -n +2 $target_proteins | cut -d ',' -f1 | uniq)

# get data for each domain
echo "get data for each domain"

echo "pfamA_acc;pfamA_id;clan_id;clan_acc;comment;selected" > domain_info.csv

i=0
for domain in $domains
do  
    echo $domain
    line=$( mysql -h $host -u $user -p$password -P $port pfam_35_0 \
    --quick -e "SELECT pfamA.pfamA_acc,pfamA.pfamA_id,clan.clan_id,clan.clan_acc,pfamA.comment FROM clan,clan_membership,pfamA WHERE (clan.clan_acc = clan_membership.clan_acc AND clan_membership.pfamA_acc = pfamA.pfamA_acc AND pfamA.pfamA_acc = '$domain') OR (pfamA.pfamA_acc = '$domain')" \
    | tail -n 1 | tr '\t' ';')
    if [ $i -lt 50 ];
        then
        line=$line";True"
    fi

    echo $line
    echo $line >> data/domain_info.csv
    i=$((i+1))
  
done
