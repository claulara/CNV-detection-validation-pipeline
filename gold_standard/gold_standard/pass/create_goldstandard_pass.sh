#!/bin/bash

source /opt/mambaforge/mambaforge/etc/profile.d/conda.sh
conda activate cnvs



###################### STUDY PREPROCESSING ######################

# Extract the first four columns and sort the variants for each study
inputs=("/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/benchmarks/1000genomes/PASS/1000genomes_pass_filtered.bed"\
        "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/benchmarks/giab/PASS/giab_pass_filtered.bed"\
        "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/benchmarks/metasv/PASS/metasv_pass_filtered.bed"\
        "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/benchmarks/svclassify/svclassify_filtered.bed")

studies_pass_dir="/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/pass/estudios/"
mkdir -p "$studies_pass_dir"

# Process each study individually
for input_file in "${inputs[@]}"; do
    filename=$(basename "$input_file")
    cut -f1-4 "$input_file" | sort -k1,1V -k2,2n > "$studies_pass_dir/$filename"
    echo "Columns (chrom start end svtype) extracted and sorted individually for overlapping variants in study $filename"
done



###################### LAX MERGE FUSION ######################

# Create output directory for the fusion process
fusion_dir="/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/pass/fusion_merge/"
mkdir -p "$fusion_dir"

# Combine and sort BED files
# At this stage, all studies must be sorted in the same format
combined_sorted_bed="${fusion_dir}/estudios_pass_sorted.bed"
cat "${studies_pass_dir}"/*.bed | bedtools sort > "$combined_sorted_bed"

# Merge with variant collapse to include SV type
merged_bed="${fusion_dir}/fusion_merge.bed"
bedtools merge -i "$combined_sorted_bed" -c 4 -o collapse > "$merged_bed"
awk -F'\t' '{if (gsub(",", "&", $4) >= 1) print $0}' OFS='\t' "${merged_bed}" > "${fusion_dir}/fusion_merge_2ormore.bed"
awk -F'\t' '{if (gsub(",", "&", $4) >= 2) print $0}' OFS='\t' "${merged_bed}" > "${fusion_dir}/fusion_merge_3ormore.bed"
awk -F'\t' '{if (gsub(",", "&", $4) >= 3) print $0}' OFS='\t' "${merged_bed}" > "${fusion_dir}/fusion_merge_4.bed"


###################### STRICT MULTIINTERSECT FUSION ######################

# Create output directory for the fusion process
fusion_dir="/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/pass/fusion_intersect/"
mkdir -p "$fusion_dir"

# Multiintersect 
intersect_bed="${fusion_dir}/fusion_intersect.bed"
bedtools  multiinter -i "${studies_pass_dir}"/*.bed -f 0.5 -r > "$intersect_bed"
awk '$4 >= 2' "${fusion_dir}/fusion_intersect.bed" | bedtools merge -i - > "${fusion_dir}/fusion_intersect_2ormore.bed" 
awk '$4 >= 3' "${fusion_dir}/fusion_intersect.bed" | bedtools merge -i - > "${fusion_dir}/fusion_intersect_3ormore.bed"
awk '$4 >= 4' "${fusion_dir}/fusion_intersect.bed" | bedtools merge -i - > "${fusion_dir}/fusion_intersect_4.bed"
