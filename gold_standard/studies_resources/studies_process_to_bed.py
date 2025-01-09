import gzip
import re

def open_file(file_path, mode='rt'):
    """
    Opens a file for reading or writing, depending on whether it is compressed (.gz) or not.
    """
    if file_path.endswith(".gz"):
        return gzip.open(file_path, mode)
    else:
        return open(file_path, mode)


def extract_sample_vcf(vcf_path, sample_id, output_vcf_path):
    """
    Extracts a specific sample from a VCF file and writes it to a new VCF file.
    """
    try:
        with open_file(vcf_path, 'rt') as input_file:
            header = None
            for line in input_file:
                if line.startswith("#CHROM"):
                    header = line.strip().split("\t")
                    break
        # check if the sample ID is in header
        if not header or sample_id not in header:
            print(f"The {sample_id} sample was not found in the VCF file.")
            return None
        column = header.index(sample_id)
        print(f"The {sample_id} sample is located in the {column} column.")
        index = list(range(8)) + [column]

        # create a new VCF file just with the specific sample
        with open_file(vcf_path, 'rt') as input_file:
            with open_file(output_vcf_path, 'wt') as output_file:
                for line in input_file:

                    # metadata lines
                    if line.startswith("##"):
                        output_file.write(line)

                    # header line
                    elif line.startswith("#CHROM"):
                        header = line.strip().split("\t")
                        new_header = [header[i] for i in index]
                        output_file.write("\t".join(new_header) + "\n")

                    # variants lines
                    else:
                        line_split = line.strip().split("\t")
                        new_line = [line_split[i] for i in index]
                        # ignore those variants with 0|0 genotype (homozygous reference)
                        if new_line[-1] != "0|0":
                            output_file.write("\t".join(new_line) + "\n")

        return output_vcf_path

    except FileNotFoundError:
        print(f"Sample extraction not performed because file '{vcf_path}' was not found. Please check the path.")
        return None


def process_vcf_quality(vcf_path, pass_vcf_path, lowquality_vcf_path):
    """
    Split a VCF file into two separate files based on the FILTER field.
    Variants with 'PASS' (high quality) in the FILTER field go to one file, and others (low quality) go to another.
    """
    
    main_chrom = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
    try:
        with open_file(vcf_path, 'rt') as input_file:
            with open_file(pass_vcf_path, 'wt') as pass_file, open_file(lowquality_vcf_path, 'wt') as lowquality_file:
                header_written = False
                filter_column_index = None

                for line in input_file:

                    # metadata lines, same for both new VCF
                    if line.startswith("##"):
                        pass_file.write(line)
                        lowquality_file.write(line)

                    # header line, same for both new VCF
                    elif line.startswith("#CHROM"):
                        pass_file.write(line)
                        lowquality_file.write(line)
                        header = line.strip().split("\t")
                        # find the index for 'FILTER' column
                        filter_column_index = header.index("FILTER")
                        header_written = True

                    # variants lines    
                    else:
                        if filter_column_index is not None:
                            line_split = line.strip().split("\t")
                            filter_value = line_split[filter_column_index]

                            # write a different VCF depending on the 'FILTER' field
                            if filter_value == "PASS":
                                pass_file.write("\t".join(line_split) + "\n")
                            else:
                                lowquality_file.write("\t".join(line_split) + "\n")
                        else:
                            print("FILTER column not found in header.")

    except FileNotFoundError:
        print(f"Subset creeation not performed because file '{vcf_path}' was not found. Please check the path.")


def vcf_to_bed(vcf_path, bed_path):
    """
    Converts a VCF file to a BED file by extracting relevant information.

    """
    main_chrom = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
    try:

        with open_file(vcf_path, 'rt') as infile, open(bed_path, 'wt') as outfile:
            # optionally write a header to the BED file
            # outfile.write("chrom\tpos\tend\tsvtype\tlength\tqual\tcallmethod\tgenotype\n")
            na12878_index = 0
            for line in infile:
                # skip metadata and header lines
                if line.startswith('##'):
                    continue
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    if 'NA12878' in header:
                        na12878_index = header.index('NA12878')
                    continue

                # columns and information to keep
                columns = line.strip().split('\t')
                chrom = columns[0]

                # add 'chr' prefix if missing
                if "chr" not in chrom: 
                    chrom = f"chr{chrom}"

                # remove alternatie chroms from GIAB study
                if chrom in main_chrom: 
                    # convert pos from 1-based to 0-based for beds
                    pos = int(columns[1]) - 1  
                    qual = columns[6]
                    info = columns[7]
                    genotype = columns[na12878_index]

                    # extract 'END' position from INFO field
                    try:
                        end_match = re.search(r'\bEND=(\d+)', info)
                        end = end_match.group(1) if end_match else ''
                        if end in {'0', ''}:
                            end_match = re.search(r'\bSVLEN=(\d+)', info) 
                            end = (int(end_match.group(1)) + int(pos)) if end_match else ''
                    except IndexError:
                        continue

                    # SV length
                    length = int(end) - int(pos)

                    # Extract the variant type 'SVTYPE' from the INFO field
                    try:
                        svtype = info.split('SVTYPE=')[1].split(';')[0]
                    except IndexError:
                        svtype = ''

                    # extract call method from the INFO field
                    if 'SVMETHOD=' in info:
                        callmethod = info.split(';SVMETHOD=')[1].split(';')[0]
                    elif 'CS=' in info:
                        callmethod = info.split('CS=')[1].split(';')[0]
                    elif 'NS=' in info:
                        callmethod = info.split('NS=')[1].split(';')[0]


                    outfile.write(f"{chrom}\t{pos}\t{end}\t{svtype}\t{length}\t{qual}\t{callmethod}\t{genotype}\n")

    except FileNotFoundError:
        print(f"Bed not created because the file '{vcf_path}' was not found. Please check the path.")



######################################################################################################################################################################################################

# Exctract NA12878 sample from Phase3 study 
milgenomes_filteredsample = extract_sample_vcf(vcf_path = "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standar/HG001/benchmarks/1000genomes/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.modified.vcf.gz", 
                                            sample_id = "NA12878",
                                            output_vcf_path = "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standar/HG001/benchmarks/1000genomes/NA12878_1000genomes.vcf.gz")
vcf_paths = {
    "giab": "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standar/HG001/benchmarks/giab/NA12878.sorted.vcf.gz",
    "metasv": "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standar/HG001/benchmarks/metasv/NA12878_svs.vcf.gz",
    "1000genomes": milgenomes_filteredsample
}

for study, vcf in vcf_paths.items():

    # Create beds
    vcf_to_bed(vcf, f"/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standar/HG001/benchmarks/{study}/{study}.bed")

    # Generate pass and lowqual vcfs bed of initial vcfs
    pass_vcf = f"/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standar/HG001/benchmarks/{study}/PASS/{study}_pass.vcf.gz"
    lowquality_vcf = f"/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standar/HG001/benchmarks/{study}/LowQual/{study}_lowqual.vcf.gz"
    process_vcf_quality(vcf, pass_vcf, lowquality_vcf)

    # Create pass and lowqual beds
    bed_pass = f"/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standar/HG001/benchmarks/{study}/PASS/{study}_pass.bed"
    vcf_to_bed(pass_vcf, bed_pass)
    bed_lowquality = f"/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standar/HG001/benchmarks/{study}/LowQual/{study}_lowqual.bed"
    vcf_to_bed(lowquality_vcf, bed_lowquality)


######################################################################################################################################################################################################
