# Making a tree using iqtree  

First, we need to convert snps.vcf to fasta - vcf2phylip.py python program

```
wget https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py
chmod +x vcf2phylip.py

module add python

python3 vcf2phylip.py -i /storage/brno12-cerit/home/alena_bartonova/RAD_Paph_Lim/camilla_gstacks/populations_output_geno20_maf05_all/populations.snps.vcf -f
```
- it makes a fasta file called populations.snps.min4.fasta

make a bash script for iqtree

```
#!/bin/bash
#PBS -N iqtree_modeltest
#PBS -l select=1:ncpus=8:mem=16gb:scratch_local=5gb
#PBS -l walltime=08:00:00

module add iqtree

cd "$SCRATCHDIR" || exit 1
cp /storage/brno2/home/alena_bartonova/populations.snps.min4.fasta .


iqtree -s populations.snps.min4.fasta -st DNA -m MFP -bb 1000 -alrt 1000 -nt AUTO -ntmax 8

echo "Output files:"
ls -lh populations.snps.min4.fasta.*
cp populations.snps.min4.fasta.* /storage/brno2/home/alena_bartonova/
```

we need as complete snp table as possible (not too much missing data)
