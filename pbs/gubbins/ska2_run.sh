#PBS -l select=1:ncpus=16:mem=32gb
#PBS -l walltime=48:00:00

## Script to run SKA 2 on the Imperial HPC 

eval "$(/rds/general/user/jd2117/home/miniforge3/bin/conda shell.bash hook)"

source /rds/general/user/jd2117/home/miniforge3/bin/activate panaroo_env


if [ $# -ne 3 ]
then 
    echo "Need 3 arguments"
    echo "Script to run ska2 via parallel on large collections of isolates"
    echo ""
    echo "Usage: "
    echo "qsub qsub ska2_run.pbs <list_of_fasta_files> <full_path_to_reference> <output_dir>"
    echo ""
    exit
else
    FILE_LIST=$1
    export REFERENCE_LOC=$2 
fi

mkdir -p $3
cd $3

START_TIME=$SECONDS
cat ${FILE_LIST} | parallel -j 16 "NGSID=\$(basename {} | sed 's/_g.*fna/_temp_skf/g');
 ska build -o \$NGSID -k 31 {} ;
 ska weed --ambig-mask \${NGSID}.skf;
MAPOUT=\$(basename {} | sed 's/_g.*fna/_aln/g');
 ska map $REFERENCE_LOC \${NGSID}.skf -o \${MAPOUT}.aln"

END=$(( SECONDS - START_TIME ))
echo "Finished in ${END} (s) "

