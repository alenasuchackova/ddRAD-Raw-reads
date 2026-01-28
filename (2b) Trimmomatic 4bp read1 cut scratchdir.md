# Trimmomatic 4bp cut from read1 using scratchdir

```
#PBS -N trimmomatic_4bp_loop
#PBS -l select=1:ncpus=3:mem=70gb:scratch_local=70gb
#PBS -l walltime=05:00:00

PREFIX="Apat2"
RAWDIR="/storage/brno12-cerit/home/alena_bartonova/RAD3_Apat2_red_sap_riv/trimmomatic_adaptercut"
OUTDIR="/storage/brno12-cerit/home/alena_bartonova/RAD3_Apat2_red_sap_riv/trimmomatic_4bpcut"

trap 'clean_scratch' TERM EXIT
cd "$SCRATCHDIR" || exit 1


cp ${RAWDIR}/${PREFIX}*_1.paired.fq.gz "$SCRATCHDIR"/


module add trimmomatic


for pool in ${PREFIX}*_1.paired.fq.gz; do
    echo "Processing $pool"
    
    pool_name=${pool%_1.paired.fq.gz}
    output_file="${pool_name}_1.paired_4bpcut.fq.gz"

    trimmomatic SE -threads 3 "$pool" "$output_file" HEADCROP:4
done

cp *_1.paired_4bpcut.fq.gz "$OUTDIR"/
```

Note: 3 cpu and 70gb mem: memory used around 50 gb, 39% of cpu time, took 2 h 30 min
