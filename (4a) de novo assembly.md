# de novo assembly using stacks  

- use instead of steps 4-6 when a reference genome is not available  
-  we need demultiplexed samples (fastq) and a population map

```
#PBS -N de_novo_iris
#PBS -l select=1:ncpus=10:mem=100gb:scratch_local=150gb
#PBS -l walltime=20:00:00

trap 'clean_scratch' TERM EXIT
cd "$SCRATCHDIR" || exit 1

cp -r /storage/brno12-cerit/home/alena_bartonova/RAD2_Paph_Apat/demultiplex_apatura1_corr ./apatura_denovo
cp /storage/brno12-cerit/home/alena_bartonova/RAD2_Paph_Apat/demultiplex_apatura1_corr/popmap_apat1.txt .
mkdir -p denovo_output

module add stacks


denovo_map.pl -T 10 --samples apatura_denovo \
--popmap popmap_apat1.txt \
-o denovo_output \
--paired \
-X "populations: --structure --plink --vcf" -X "ustacks: -M 3 --force-diff-len"

cp -r denovo_output /storage/brno12-cerit/home/alena_bartonova/RAD2_Paph_Apat/apatura1_denovo
```
Note: 95 samples took ~10 hrs, memory used well, but only 20% cpu
