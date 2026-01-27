# Process_radtags script generator
This script is generating a set of executable scripts, which can be submitted to Metacentrum, each demultiplexing a different pool  
- changing library, adjust the PREFIX and path to barcodes and OUTDIR 
- generates a directory called pool_scripts
  
after generating scripts use a for loop to qsub them:  

` for script in pool_scripts/*.sh; do qsub "$script"; done `  

Note: one pool took 30-60 mins, 3 cpu used well, but only about 25% of memory (50gb)

```
#!/bin/bash

# Output directory for generated scripts
GEN_DIR="./pool_scripts"
mkdir -p "$GEN_DIR"

# Variables
PREFIX="Apat2"
RAWDIR="/storage/brno12-cerit/home/alena_bartonova/RAD3_Apat2_red_sap_riv"
TRIM_4BP_DIR="${RAWDIR}/trimmomatic_4bpcut"
TRIM_ADAPTER_DIR="${RAWDIR}/trimmomatic_adaptercut"
BARCODEDIR="${RAWDIR}/demultiplex_apatura2/barcodes_Apat2"
OUTDIR="${RAWDIR}/demultiplex_apatura2" 

POOL_NUMBERS=(011 179 220 256 347 455 509 526 604 625)

for POOL in "${POOL_NUMBERS[@]}"; do
    SCRIPT="${GEN_DIR}/run_${PREFIX}_SP${POOL}.sh"

    cat << EOF > "$SCRIPT"
#!/bin/bash
#PBS -N ${PREFIX}_SP${POOL}
#PBS -l select=1:ncpus=3:mem=50gb:scratch_local=50gb
#PBS -l walltime=5:00:00
#PBS -j oe

module add stacks

trap 'clean_scratch' TERM EXIT
cd "\$SCRATCHDIR" || exit 1

echo "Processing pool ${PREFIX}_SP${POOL}"

# Copy input files
cp "${TRIM_4BP_DIR}/${PREFIX}_SP${POOL}_1.paired_4bpcut.fq.gz" .
cp "${TRIM_ADAPTER_DIR}/${PREFIX}_SP${POOL}_2.paired.fq.gz" .
cp "${BARCODEDIR}/barcodes_SP${POOL}.txt" .

# Run process_radtags
process_radtags -i gzfastq \\
    -1 "${PREFIX}_SP${POOL}_1.paired_4bpcut.fq.gz" \\
    -2 "${PREFIX}_SP${POOL}_2.paired.fq.gz" \\
    -b "barcodes_SP${POOL}.txt" \\
    -o "\$SCRATCHDIR" \\
    --renz-1 ecoRI --renz-2 mseI --threads 3 -c -q -r

cp "\$SCRATCHDIR"/*.fq.gz "$OUTDIR/"

echo "Finished pool ${PREFIX}_SP${POOL}"
EOF

    chmod +x "$SCRIPT"
    echo "Created $SCRIPT"
done

```
