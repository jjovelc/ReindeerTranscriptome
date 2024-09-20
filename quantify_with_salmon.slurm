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
SIF=salmon_singularity.sif
TRANSCRIPTOME_IDX=bt-rt_transcriptome_idx
TGMAP=${DIR}/"txp2gene.tsv"

# Construct the lists of barcode and read files
BARCODES_FILES=()
READS_FILES=()
for BARCODES in ${DIR}/*_R1.fq.gz; do
    BARCODES_FILES+=($BARCODES)
    READS="${BARCODES%_R1.fq.gz}_R2.fq.gz"
    READS_FILES+=($READS)
done

# Run Salmon Alevin for each pair
for i in "${!BARCODES_FILES[@]}"; do
    OUTPUT_DIR="${DIR}/alevin_output_240918$(basename ${BARCODES_FILES[$i]} _R1.fq.gz)"
    
    singularity exec "${DIR}/$SIF" salmon alevin -l ISR \
        -i $DIR/$TRANSCRIPTOME_IDX \
        -1 ${BARCODES_FILES[$i]} \
        -2 ${READS_FILES[$i]} \
        -o $OUTPUT_DIR \
        --tgMap $TGMAP \
        --dumpFeatures \
        --dumpUmiGraph \
        --chromium \
        --dumpMtx \
        -p 10
done

end_time=$(date +%s)
echo "Job ended on: $(date)"

# Calculate the difference in times to get the duration
duration=$((end_time - start_time))

hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))
echo "Job duration: $hours hours, $minutes minutes, $seconds seconds"
