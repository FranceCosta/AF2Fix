#! /bin/bash


input=output/proteins_per_domain.csv
output=output/target_proteins.tsv
number_of_proteins=10

while getopts ":i:o:n:u:p:h:P:" flag; do
    case $flag in
        i) input=${OPTARG};;
        o) output=${OPTARG};;
        n) number_of_proteins=${OPTARG};;
        u) SQLuser=${OPTARG};; # user for pfam MySQL
        p) SQLpassword=${OPTARG};; # password for pfam MySQL
        h) SQLhost=${OPTARG};; # host for pfam MySQL
        P) SQLport=${OPTARG};; # port for pfam MySQL
    esac
done

echo pfamA_acc uniprot_acc seq_start seq_end model_start model_end domain_bits_score domain_evalue_score sequence_bits_score sequence_evalue_score sequence > $output
while read row
do
  echo 'row: '$row
  arr=(${row//,/ })
  pfamA_acc=${arr[0]}
  successfull_proteins=0
  iteration=0
  for pfamseq_acc in ${arr[@]:1}
  do
    #pfamseq_acc=$(echo ${pfamseq_acc//_/ } | awk '{print $1}')
    echo 'pfamseq_acc: '$pfamseq_acc' success: '$successfull_proteins' iteration: '$iteration
    result=$(mysql -h $SQLhost -u $SQLuser -p$SQLpassword -P $SQLport pfam_35_0 --quick \
    -e "SELECT pfamA_acc,uniprot_reg_full.uniprot_acc,seq_start,seq_end,model_start,model_end,domain_bits_score,domain_evalue_score,sequence_bits_score,sequence_evalue_score,uniprot.sequence\
    from uniprot_reg_full,uniprot WHERE pfamA_acc='$pfamA_acc' AND (uniprot_reg_full.uniprot_acc='$pfamseq_acc' OR uniprot.uniprot_id='$pfamseq_acc') AND uniprot_reg_full.uniprot_acc = uniprot.uniprot_acc" \
    | tail -n +2)
    echo $result
    if [[ ! -z "$result" ]]
    then
      echo $result >> $output
      successfull_proteins=$(( successfull_proteins + 1 ))
    fi
      iteration=$((iteration + 1))

    if [ $successfull_proteins -ge $number_of_proteins ]
    then
      echo 'breaking'
      break 1
    fi
    done

    if ((successfull_proteins < $number_of_proteins))
    then
      echo 'pfamA: '$pfamA_acc' success:' $successfull_proteins
    fi
done < $input
