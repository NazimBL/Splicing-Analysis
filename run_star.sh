cd ~/macaca_ref/GSE_test

for sample in SRR19357307 SRR19357308 SRR19357309 SRR19357310
do
  STAR --runThreadN 8 \
       --genomeDir ~/macaca_ref/star_index \
       --readFilesIn ${sample}_1.fastq ${sample}_2.fastq \
       --outFileNamePrefix ${sample}_ \
       --outSAMtype BAM SortedByCoordinate
done
