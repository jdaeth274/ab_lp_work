#PBS -l select=1:ncpus=16:mem=32gb
#PBS -l walltime=48:00:00
#PBS -J 1-6
## Script to run SKA 2 on the Imperial HPC 

eval "$(/rds/general/user/jd2117/home/miniforge3/bin/conda shell.bash hook)"

source /rds/general/user/jd2117/home/miniforge3/bin/activate ska2_env

CURRENT_FILE=$(head -n $PBS_ARRAY_INDEX /rds/general/user/jd2117/home/acba_legion_2024/ska_input/acba_files_references.txt | tail -n 1)
FILE_LIST=$(echo $CURRENT_FILE | awk '{print $1}')
REFERENCE_LOC=$(echo $CURRENT_FILE | awk '{print $2}')

export $REFERENCE_LOC 

OUTPUT_NAME=$(basename $FILE_LIST | sed 's/\..*$/_ska_res/g')

mkdir -p "~/../ephemeral/acba_legion/ska2_res/${OUTPUT_NAME}/"
cd "~/../ephemeral/acba_legion/ska2_res/${OUTPUT_NAME}/"

START_TIME=$SECONDS
cat ${FILE_LIST} | parallel -j 16 "NGSID=\$(basename {} | sed 's/_g.*fna/_temp_skf/g');
 ska build -o \$NGSID -k 31 {} ;
 ska weed --ambig-mask \${NGSID}.skf;
MAPOUT=\$(basename {} | sed 's/_g.*fna/_aln/g');
 ska map $REFERENCE_LOC \${NGSID}.skf -o \${MAPOUT}.aln"

END=$(( SECONDS - START_TIME ))
echo "Finished in ${END} (s) "

