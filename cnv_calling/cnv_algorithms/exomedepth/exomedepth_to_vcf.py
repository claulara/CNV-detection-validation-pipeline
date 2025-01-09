#!/usr/bin/python3.8

import optparse
import re, os, csv, gzip 

def normalize_variants(variant):
    """
    Normalize variant names
    """
    if variant in ['duplication', 'dup', 'DUP', '3', '4', '5', 'CN3', 'CN4', 'CN5', 'CN6', 'CN7', 'CN8', 'CNV gain', 'gain']:
        return 'DUP' 
    elif variant in ['deletion', 'del', 'DEL', 'CNV loss', 'loss', '0', 'CN0', '1', 'CN1']:
        return 'DEL'  
    elif variant in ['insertion', 'ins', 'INS']:
        return 'INS' 
    elif variant in ['inversion', 'inv', 'INV']:
        return 'INV'  
    else:
        return 'CNV'  


def bed_to_vcf(input_bed, output_vcf):
    """
    Converts a BED file with structural variant information into a VCF file.
    """
    sample_name = os.path.basename(input_bed).split('.')[0]
    
    vcf_header = f"""##fileformat=VCFv4.2
##fileformat=VCFv4.1
##fileDate=20241001
##source=GenerateSVCandidates 1.6.0
##reference=file:///ingemm/ref/Ref/hg19/fasta/ucsc.hg19.fasta
##contig=<ID=chrM,length=16571>
##contig=<ID=chr1,length=249250621>
##contig=<ID=chr2,length=243199373>
##contig=<ID=chr3,length=198022430>
##contig=<ID=chr4,length=191154276>
##contig=<ID=chr5,length=180915260>
##contig=<ID=chr6,length=171115067>
##contig=<ID=chr7,length=159138663>
##contig=<ID=chr8,length=146364022>
##contig=<ID=chr9,length=141213431>
##contig=<ID=chr10,length=135534747>
##contig=<ID=chr11,length=135006516>
##contig=<ID=chr12,length=133851895>
##contig=<ID=chr13,length=115169878>
##contig=<ID=chr14,length=107349540>
##contig=<ID=chr15,length=102531392>
##contig=<ID=chr16,length=90354753>
##contig=<ID=chr17,length=81195210>
##contig=<ID=chr18,length=78077248>
##contig=<ID=chr19,length=59128983>
##contig=<ID=chr20,length=63025520>
##contig=<ID=chr21,length=48129895>
##contig=<ID=chr22,length=51304566>
##contig=<ID=chrX,length=155270560>
##contig=<ID=chrY,length=59373566>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">
##INFO=<ID=CIGAR,Number=A,Type=String,Description="CIGAR alignment for each alternate indel allele">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}"""

    with open(input_bed, "r") as bed, open(output_vcf, "w") as vcf:
        vcf.write(vcf_header.strip() + "\n")  
        input_file = csv.reader(bed, delimiter="\t")
    
        for line in input_file:
            chrom, start, end, svtypes = line[:4] 
            # add "chr" prefix if missing
            if not chrom.startswith("chr"):  
                chrom = f"chr{chrom}"

            start = int(start)
            # convert pos from 0-based to 1-based for vcfs
            pos = start + 1 
            end = int(end)
            svlen = abs(end - start)
            
            for svtype in svtypes.split(","):
                svtype = svtype.strip() 
                info = f"SVTYPE={svtype};END={end};SVLEN={svlen}"
                
                vcf.write(f"{chrom}\t{pos}\t.\tN\t<{svtype}>\t.\tPASS\t{info}\t.\t.\n")


def exomedepth_to_bed(file_path, bed_path):
    """
    Conversion of exomedepth output to normalized BED and VCF files
    """

    with open(file_path, 'rt') as infile, open(bed_path, 'wt') as outfile:

        # remove header
        next(infile)
        for line in infile:

            # columns to keep
            columns = [col.strip('"') for col in line.strip().split(',')]
            chrom = columns[7]
            start = int(columns[5]) 
            end = int(columns[6]) 
            svtype = normalize_variants(columns[3])  
            length = int(end) - int(start)

            # Escribir la línea en formato BED
            outfile.write(f"{chrom}\t{start}\t{end}\t{svtype}\t{length}\n")

    bed_to_vcf(bed_path, bed_path.replace(".bed", ".vcf"))
    
    return


def run(argv=None):

    parser = optparse.OptionParser()
    parser.add_option('--f', default=None, help='Manta output path', dest='file_path')
    parser.add_option('--o', default=None, help='Normalized manta BED path', dest='bed_path')
    (options, args) = parser.parse_args(argv[1:])

    file_path = options.file_path
    bed_path = options.bed_path
    exomedepth_to_bed(file_path, bed_path)


    print(f'BED y VCF creado')


if __name__ == '__main__':
    import sys
    run(sys.argv)
