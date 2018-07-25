#!/usr/bin/env python

import sys, pysam, argparse

parser = argparse.ArgumentParser(
        usage="%(prog)s barcodefile bamfile",
        description="Filter bamfile to keep alignment with a valid cell barcode")
parser.add_argument( "bamfile", type=str,
            help="Path to the SAM/BAM files, must be indexed")
parser.add_argument( "barcodefile", type=str,
            help="barcodefile, output from cellranger")

args = parser.parse_args()
barcodefile = args.bamfile
baminputfile = args.barcodefile

barcodes = []
counter = 0

for line in open(barcodefile):
    barcodes.append(line.strip())

bamfile = pysam.AlignmentFile(baminputfile, "rb")
fnm = "filter."+sys.argv[2]
outfile = pysam.AlignmentFile(fnm, "wb", template=bamfile)
for read in bamfile.fetch():
    try:
        bc = read.get_tag('CB')
        if bc in barcodes:
            outfile.write(read)
    except: continue

outfile.close()
bamfile.close()
