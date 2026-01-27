# Trimmomatic adaptercut scratch
uses scratchdir on Metacentrum and can run only part of the data by changing PREFIX

```
#PBS -N trimmomatic_adapterclip_loop
#PBS -l select=1:ncpus=5:mem=100gb:scratch_local=100gb
#PBS -l walltime=10:00:00

module add trimmomatic

PREFIX="Apat2"
RAWDIR="/storage/brno12-cerit/home/alena_bartonova/RAD3_Apat2_red_sap_riv/trimmomatic_adaptercut/raw"
OUTDIR="/storage/brno12-cerit/home/alena_bartonova/RAD3_Apat2_red_sap_riv/trimmomatic_adaptercut"
ADAPTERS="/storage/brno12-cerit/home/alena_bartonova/RAD3_Apat2_red_sap_riv/trimmomatic_adaptercut/raw/adapters.fas"

trap 'clean_scratch' TERM EXIT
cd "$SCRATCHDIR" || exit 1



cp ${RAWDIR}/${PREFIX}*_1.fq.gz "$SCRATCHDIR"/
cp ${RAWDIR}/${PREFIX}*_2.fq.gz "$SCRATCHDIR"/
cp "$ADAPTERS" "$SCRATCHDIR"/    

for R1 in ${PREFIX}*_1.fq.gz; do
    R2="${R1%_1.fq.gz}_2.fq.gz"

    trimmed_R1="$R1"
    trimmed_R2="$R2"

    trimmed_R1="${trimmed_R1%.*.*}"
    trimmed_R2="${trimmed_R2%.*.*}"

    out_R1_paired="${trimmed_R1%.*}.paired.fq.gz"
    out_R1_unpaired="${trimmed_R1%.*}.unpaired.fq.gz"
    out_R2_paired="${trimmed_R2%.*}.paired.fq.gz"
    out_R2_unpaired="${trimmed_R2%.*}.unpaired.fq.gz"

    trimmomatic PE -threads 5 "$R1" "$R2" \
        "$out_R1_paired" "$out_R1_unpaired" \
        "$out_R2_paired" "$out_R2_unpaired" \
        ILLUMINACLIP:${SCRATCHDIR}/adapters.fas:2:30:10 \
        MINLEN:80
done


cp *.paired.fq.gz "$OUTDIR"/
cp *.unpaired.fq.gz "$OUTDIR"/

```
note: with 4 threads 100 gb memory, 10 pools took cca 2 h 30 min, 65% cpu usage - better use less next time
