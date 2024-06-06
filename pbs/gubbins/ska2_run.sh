#PBS -l select=1:ncpus=1:mem=16gb
#PBS -l walltime=24:00:00
#PBS -J 1-10
## Script to run SKA 2 on the Imperial HPC for the acba gc2 lineage

source /rds/general/user/jd2117/home/miniforge3/bin/activate ska2_env

CURRENT_FILE=$(head -n $PBS_ARRAY_INDEX /rds/general/user/jd2117/home/acba_legion_2024/ska_input/gc2_splits.txt | tail -n 1)

REFERENCE_LOC="/rds/general/user/jd2117/home/acba_legion_2024/ska_input/GCA_001573125.1_ASM157312v1_genomic.fna"

OUTPUT_NAME="GC2_SKA2_RES"

mkdir -p "~/../ephemeral/acba_legion/ska2_runs/gc2_ska/"
cd "~/../ephemeral/acba_legion/ska2_runs/gc2_ska/"

START_TIME=$SECONDS
cat ${CURRENT_FILE} | while read line;
do
NGSID=$(basename $line | sed 's/_g.*fna/_temp_skf/g')
ska build -o $NGSID -k 31 $line ;
ska weed --ambig-mask ${NGSID}.skf;
MAPOUT=$(basename $line | sed 's/_g.*fna/_aln/g');
ska map $REFERENCE_LOC ${NGSID}.skf -o ${MAPOUT}.aln

done
END=$(( SECONDS - START_TIME ))
echo "Finished in ${END} (s) "

