#PBS -l select=1:ncpus=16:mem=32gb
#PBS -l walltime=48:00:00

## Script to run SKA 2 on the Imperial HPC 

eval "$(/rds/general/user/jd2117/home/miniforge3/bin/conda shell.bash hook)"

source /rds/general/user/jd2117/home/miniforge3/bin/activate panaroo_env


FILE_LIST="~/acba_legion_2024/ska_input/gc2_5092_paths.txt"
export REFERENCE_LOC="~/acba_legion_2024/ska_input/GCA_001573125.1_ASM157312v1_genomic.fna" 

mkdir -p ~/../ephemeral/acba_legion/ska2_res/GC2/
cd ~/../ephemeral/acba_legion/ska2_res/GC2/

START_TIME=$SECONDS
cat ${FILE_LIST} | parallel -j 16 "NGSID=\$(basename {} | sed 's/_g.*fna/_temp_skf/g');
 ska build -o \$NGSID -k 31 {} ;
 ska weed --ambig-mask \${NGSID}.skf;
MAPOUT=\$(basename {} | sed 's/_g.*fna/_aln/g');
 ska map $REFERENCE_LOC \${NGSID}.skf -o \${MAPOUT}.aln"

END=$(( SECONDS - START_TIME ))
echo "Finished in ${END} (s) "

