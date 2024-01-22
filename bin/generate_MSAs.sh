
input_tsv=output/target_proteins.tsv
database=/nfs/research/agb/research/francesco/software/localcolabfold/databases/

while getopts ":i:d:" flag; do
    case $flag in
        i) input_tsv=${OPTARG};;
        d) database=${OPTARG};;
    esac
done

#mkdir tmp

# get sequences from file
tail -n +2 $input_tsv | cut -d ' ' -f 1,2,11 | awk -F ' ' '{print ">"$1"_"$2"\n"$3""}' > sequences.fa

# generate MSAs
#conda activate /nfs/research/agb/research/francesco/software/localcolabfold/colabfold-conda
colabfold_search --use-templates 0 --threads 96 sequences.fa $database ./

# split MSAs into domain folders
for MSA in $( ls *.a3m )
do
  domain=$( head $MSA -n 1 | tr ">" " " | tr "_" " " | awk '{ print $1 }' )
  protein=$( head $MSA -n 1 | tr ">" " " | tr "_" " " | awk '{ print $2 }' )
  mkdir -p $domain
  cp $MSA $domain/$protein.a3m
  #touch tmp/$domain
done
