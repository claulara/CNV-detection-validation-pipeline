import os
import csv

def gold_standard_to_vcf(input_bed, output_vcf):
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


GS = ["/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/pass/gold_standard_pass.bed",
    "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/pass/gold_standard_pass_IDTv1.bed",
    "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/pass/gold_standard_pass_IDTv2.bed",
    "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/pass/gold_standard_pass_ROCHEv1.bed",
    "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/general/gold_standard_general.bed",
    "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/general/gold_standard_general_IDTv1.bed",
    "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/general/gold_standard_general_IDTv2.bed",
    "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/general/gold_standard_general_ROCHEv1.bed"
    ]
for gs in GS:
    gold_standard_to_vcf(gs, gs.replace(".bed", ".vcf"))
