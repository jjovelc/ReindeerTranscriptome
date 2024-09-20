#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=48:00:00
#SBATCH --mem=8G

####### Set environment variables ###############

start_time=$(date +%s)
echo "Job started on: $(date)"

####### Run your script #########################
# Run the Salmon index command inside the Singularity container
DIR=${PWD}
TRANSCRIPTOME=bt-rt_transcriptome.fa
SIF=salmon_singularity.sif

if [[ $TRANSCRIPTOME == *.fa ]]; then
    INDEXFILE="${TRANSCRIPTOME%.fa}_idx"
elif [[ $TRANSCRIPTOME == *.fasta ]]; then
    INDEXFILE="${TRANSCRIPTOME%.fasta}_idx"
elif [[ $TRANSCRIPTOME == *.fna ]]; then
    INDEXFILE="${TRANSCRIPTOME%.fna}_idx"
fi


singularity exec "${DIR}/$SIF"  salmon index -t "${DIR}/$TRANSCRIPTOME"  -i "$INDEXFILE"

end_time=$(date +%s)
echo "Job ended on: $(date)"

# Calculate the difference in times to get the duration
duration=$((end_time - start_time))

hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))
echo "Job duration: $hours hours, $minutes minutes, $seconds seconds"
