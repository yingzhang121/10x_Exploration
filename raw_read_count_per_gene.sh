#!/bin/bash -l

# Require: python3 bamtools samtools R
# Make sure the input files are available

CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
barcodefile=barcodes.tsv
bamfile=possorted_genome_bam.bam
gfffile=genes.gtf

bamtools split -in ${bamfile} -reference

for chr in ${CHROMOSOMES}; do
    awk -v foo=${chr} '$1==foo&&$3=="gene"' ${gfffile} > genes.${chr}.gtf
    python filter_barcode.py ${barcodefile} *.REF_${chr}.bam
    python count_raw_read.py filter.*.REF_${chr}.bam ${barcodefile} genes.${chr}.gtf > ${chr}.count
done

R --vallina < count.R
