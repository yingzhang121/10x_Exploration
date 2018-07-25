#!/usr/bin/env python

# HTSeq, count.py

import sys, argparse, itertools, warnings, traceback
import HTSeq

def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2

def prepare_features( gff_filename, feature_type = "gene", id_attribute = "gene_id", 
                    quiet = "TRUE", stranded = "yes" ):
    
    features = HTSeq.GenomicArrayOfSets("auto", stranded != "no")
    gff = HTSeq.GFF_Reader(gff_filename)
    counts = {}
    i = 0
    try:
        for f in gff:
            if f.type == feature_type:
                try:
                    feature_id = f.attr[id_attribute]
                except KeyError:
                    raise ValueError("Feature %s does not contain a '%s' attribute" %
                                     (f.name, id_attribute))
                if stranded != "no" and f.iv.strand == ".":
                    raise ValueError("Feature %s at %s does not have strand information but you are "
                                     "running htseq-count in stranded mode. Use '--stranded=no'." %
                                     (f.name, f.iv))
                features[f.iv] += feature_id
                counts[f.attr[id_attribute]] = 0
            i += 1
            if i % 100000 == 0 and not quiet:
                print("%d GFF lines processed.\n" % i, file=sys.stderr)
    except:
        print(
            "Error occured when processing GFF file (%s):\n" %
            gff.get_line_number_string(), file=sys.stderr)
        raise

    if not quiet:
        print("%d GFF lines processed.\n" % i, file=sys.stderr)
    
    if len(counts) == 0:
        print("Warning: No features of type '%s' found.\n" % feature_type, file=sys.stderr)
        sys.exit("GFF file is corrupt.")

    return features, counts

def count_reads_in_features(bamfile, features, counts_all, 
                        multimapped_mode = "none", stranded = "yes",
                        samtype = "bam", overlap_mode = "union", minaqual = 10,
                        secondary_alignment_mode = "ignore",
                        supplementary_alignment_mode = "ignore", quiet = "TRUE"):

    # CIGAR match characters (including alignment match, sequence match and mismatch)
    com = ('M', '=', 'X')

    if samtype == "sam":
        SAM_or_BAM_Reader = HTSeq.SAM_Reader
    elif samtype == "bam":
        SAM_or_BAM_Reader = HTSeq.BAM_Reader
    else:
        raise ValueError("Unknown input format %s specified." % samtype)

    try:
        if bamfile != "-":
            read_seq_file = SAM_or_BAM_Reader(bamfile)
            read_seq = read_seq_file
            first_read = next(iter(read_seq))
        else:
            read_seq_file = SAM_or_BAM_Reader(sys.stdin)
            read_seq_iter = iter(read_seq_file)
            first_read = next(read_seq_iter)
            read_seq = itertools.chain([first_read], read_seq_iter)
    except:
        print(
            "Error occured when reading beginning of SAM/BAM file.\n", file=sys.stderr)
        raise

    try:
        i = 0
        for r in read_seq:
            if i > 0 and i % 100000 == 0 and not quiet:
                print("%d SAM alignment record%s processed.\n" %
                    (i, "s" if not pe_mode else " pairs"), file=sys.stderr)
            i += 1
            # check r.tag['CB']
            tags = {x:y for x,y in r.optional_fields}
            try: barcode = tags['CB']
            except: 
                print(tags, file=sys.stderr)
                raise
            if not r.aligned: continue
            if ((secondary_alignment_mode == 'ignore') and
               r.not_primary_alignment):
                continue
            if ((supplementary_alignment_mode == 'ignore') and
               r.supplementary):
                continue
            try:
                if r.optional_field("NH") > 1:
                    if multimapped_mode == 'none':
                        continue
            except KeyError:
                pass
            if r.aQual < minaqual: continue
            if stranded != "reverse":
                iv_seq = (co.ref_iv for co in r.cigar if co.type in com
                          and co.size > 0)
            else:
                iv_seq = (invert_strand(co.ref_iv)
                          for co in r.cigar if (co.type in com and
                                                co.size > 0))
            try:
                if overlap_mode == "union":
                    fs = set()
                    for iv in iv_seq:
                        if iv.chrom not in features.chrom_vectors:
                            raise UnknownChrom
                        for iv2, fs2 in features[iv].steps():
                            fs = fs.union(fs2)
                elif overlap_mode in ("intersection-strict",
                                      "intersection-nonempty"):
                    fs = None
                    for iv in iv_seq:
                        if iv.chrom not in features.chrom_vectors:
                            raise UnknownChrom
                        for iv2, fs2 in features[iv].steps():
                            if ((len(fs2) > 0) or
                               (overlap_mode == "intersection-strict")):
                                if fs is None:
                                    fs = fs2.copy()
                                else:
                                    fs = fs.intersection(fs2)
                else:
                    sys.exit("Illegal overlap mode.")
                if fs is not None and len(fs) > 0:
                    if multimapped_mode == 'none':
                        if len(fs) == 1:
                            counts_all[barcode][list(fs)[0]] += 1
                    elif multimapped_mode == 'all':
                        for fsi in list(fs):
                            counts_all[barcode][fsi] += 1
                    else:
                        sys.exit("Illegal multimap mode.")
            except: pass
    except:
        print("Error occured when processing SAM input (%s):\n" %
            read_seq_file.get_line_number_string(), file=sys.stderr)
        raise

    if not quiet:
        print("%d SAM %s processed.\n" %
            (i, "alignments " if not pe_mode else "alignment pairs"), file=sys.stderr)

    return counts_all

def get_parser():
    """ Parse input """
    desc = "Count Raw Reads mapped per gene per cell per chromosome. The counting algorithm is a simplified adoption of HTSeq function. For example, I only count the read mapping at gene level (no other feature choices)."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('bamfile', type=str, help="BAM file per chromosome")
    parser.add_argument('barcodefile', type=str, help="10x valid/filtered barcodes (cellranger output")
    parser.add_argument('gfffile', type=str, help="GFF annotation file per chromosome")
    parser.add_argument('-m', '--multimapped_mode',  type=str, default="none", choices=["none", "all"], 
        help="HTSeq - multimapped mode")
    parser.add_argument('-s', "--stranded", type=str, default="yes", choices=['yes', 'no', 'reverse'], 
        help="HTSeq - stranded, reverse = yes + reverse strand interpretation")
    parser.add_argument("-o", "--overlap_mode", type=str, default="union",
            choices=["union", "intersection-strict", "intersection-nonempty"],
            help="HTSeq - overlay mode, handle reads overlapping more than one feature") 
    parser.add_argument('-q', "--minqual", type=int, default=10, help="HTSeq - threshold for mapping quality" ) 
    
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    
    bamfile = args.bamfile
    if bamfile.endswith("bam"): ftype = "bam"
    if bamfile.endswith("sam"): ftype = "sam"
    
    barcodefile = args.barcodefile
    gfffile = args.gfffile
    
    multimode = args.multimapped_mode
    strand = args.stranded
    overlay = args.overlap_mode
    minqual = args.minqual
    #gfffile = "/scratch.global/zhan2142/starr/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf" 

    # make feature files:
    features, counts = prepare_features(gff_filename = gfffile)

    counts_all = {}
    for line in open(barcodefile):
        bc = line.strip()
        counts_all[ bc ] = counts.copy()

    results = count_reads_in_features(bamfile, features, counts_all,
                            multimapped_mode = multimode, stranded = strand, 
                            samtype = ftype, overlap_mode = overlay, minaqual = minqual)

    print("cellNo\t" + "\t".join(sorted(counts.keys())))
    print_results( results )

def print_results( results ):
    for bc in sorted(results.keys()):
        count_bc = results[bc]
        print(bc, end="\t")
        for fn in sorted(count_bc.keys()):
            print(count_bc[fn], end='\t')
        print()

main()
