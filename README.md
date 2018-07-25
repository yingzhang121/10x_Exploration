# 10x Exploration

### 10x introduction:
10x technology is a method that can sequence a large number of single cells in batch.
The basics of 10x technology is to use a two-barcode system to label each molecule. 
The first barcode is called cell barcode, and the second barcode is called UMI (unique molecule identifier).

The major demultiplex algorithm is wrapped into a software package, called cellranger.
From the website:
> Cell Ranger is a set of analysis pipelines that process Chromium single-cell RNA-seq output to align reads, 
> generate gene-cell matrices and perform clustering and gene expression analysis.

**This repository is for a few scripts I wrote for the exploration of 10x data.**

#### Extract the raw counts per gene per cell
I am wondering how the UMI count compares to the raw read count. 
So I decided to generate the count table of reads mapped per gene per cell.
I started with two outputs from cellranger: the raw bam file "possorted_genome_bam.bam" and the filtered, corrected "barcode.tsv".
Initially, I planned to split the bam file into thousands of small bam files given each cell barcode (CB). 
However, I found it is almost impractical to do so.

> There are some discussion on splitting 10x data per cell, such as
> https://www.biostars.org/p/46327/
> https://github.com/pezmaster31/bamtools/issues/135
> Basically speaking the issue is there are thousands of cells in one bam file, so neither the program
> such as bamtools can keep all output files open, nor an IO routine can do it efficiently.
> (To be frank, I tried the methods mentioned in a few posts, none worked).

Then I got suggestion to set up a RAM DISC, which I didn't pursue further because I got a third recommedation,
which is to run read count while I go through the cell barcode. In other words, I set up a hashtable (like a 2-D
array, with column to be cell, and row to be gene) and update the number of mapped reads as the bam file is processed. 

I thought it is a good idea, and chose HTSeq as the python library that handles the counting part.
* I do wish HTSeq developer could re-write their count.py script, so that I could import it and use the functions as modules *

In order to test my code, I split the raw possorted_genome_bam.bam file per chromosome. Then I realized this is a good 
way for parallelization. However, when I tested my code, I found it runs very quick, which means, Aha, a simple loop can have 
my job done within one or two hours. (it seems no need for parallelization for now, but way to go here).

So I came up with the "raw_read_count_per_gene.sh" routine to generate the raw count table.

Follow-up analysis:
------
I ran some correlation analysis between my raw read count table and the cellranger UMI count table.
They correlate very well, with pearson correlation coefficient ranging from 0.8 - 0.9 per cell with or without all those "0" counts.
Of course, the p value is very significant. However, for genes with low counts, I did observe greater variability between 
raw read counts and UMI counts.  It will be very interesting to dig it out more here.
