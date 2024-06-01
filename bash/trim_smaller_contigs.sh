#!/bin/bash
## Script to remove contigs of set length from assemblies
set -e 
START=$SECONDS
if [ $# -lt 1 ]
then
    printf "\n\
Script to remove contigs of length x from a list of  \n\
assembly files supplied, outputs to current directory \n\
Requires Bio::SeqIO in perl to be installed  
Usage: \n\
bash fastaq_limit_500_bp.sh <list_of_fastas[required]> <length_cutoff[default=500]> \n\
\n\
Make sure you use this as above \n\
"
    exit 
fi 


 FASTA_LIST=$1

## Set up the out dir

if [ $# -eq 2 ]
then 
    CUTOFF=$2
else
    CUTOFF="500"
    
fi

## Getting the dir loc for the other scripts

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
parentdir="$(dirname "$DIR")"
#pythondir="${parentdir}/python/"
perldir="${parentdir}/perl/"



## Enter the loop to remove 
echo "Entering the loop to remove contigs of length <= ${CUTOFF} "
ISO_NUM="0"
TOT_ISO=$(wc -l < $FASTA_LIST)
while read line <&3
do
CURRENT_FASTA=$(basename $line)
CURRENT_EXT=$(echo "${CURRENT_FASTA##*.}")
CURRENT_TRIM=$(echo $CURRENT_FASTA | sed "s/\.${CURRENT_EXT}/_trimmed_${CUTOFF}bp_contigs\.${CURRENT_EXT}/g")

## Run the trim from the perl script in the same parent diretory 
PERLSCRIPT="${perldir}filter_contigs.pl"

touch $CURRENT_TRIM

perl $PERLSCRIPT $CUTOFF $line $CURRENT_TRIM

ISO_NUM=$(( ISO_NUM + 1 ))

printf "\rFinished on isolate %s of %s isolates" $ISO_NUM $TOT_ISO

done 3< $FASTA_LIST

END=$(( SECONDS - START ))

printf "\nTook this long to subset %s isolates: %s (S)\n" $TOT_ISO $END


