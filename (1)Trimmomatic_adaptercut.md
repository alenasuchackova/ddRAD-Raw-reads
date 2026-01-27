#PBS -N trimmomatic_adapterclip_loop
#PBS -l select=1:ncpus=3:mem=50gb:scratch_local=50gb
#PBS -l walltime=10:00:00


module add trimmomatic

for R1 in /storage/brno12-cerit/home/alena_bartonova/RAD_Paph_Lim/Trimmomatic_adaptercut/raw/raw_*_1.fq.gz; do
    R2="${R1%_1.fq.gz}_2.fq.gz"
    trimmed_R1="${R1/raw_/}"
    trimmed_R2="${R2/raw_/}"


        trimmed_R1="${trimmed_R1%.*.*}"
        trimmed_R2="${trimmed_R2%.*.*}"


    trimmomatic PE -threads 3 "$R1" "$R2" \
        "${trimmed_R1%.*}.paired.fq.gz" "${trimmed_R1%.*}.unpaired.fq.gz" \
        "${trimmed_R2%.*}.paired.fq.gz" "${trimmed_R2%.*}.unpaired.fq.gz" \
        ILLUMINACLIP:/storage/brno12-cerit/home/alena_bartonova/RAD_Paph_Lim/Trimmomatic_adaptercut/adapters.fas:2:30:10 \
        MINLEN:80

done
