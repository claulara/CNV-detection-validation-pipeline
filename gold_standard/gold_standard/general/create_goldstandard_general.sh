#!/bin/bash

source /opt/mambaforge/mambaforge/etc/profile.d/conda.sh
conda activate cnvs


###################### STUDY PREPROCESSING ######################

# Extract the first four columns and merge variants for each study
inputs=("/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/benchmarks/1000genomes/LowQual/1000genomes_lowqual_filtered.bed"\
        "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/benchmarks/giab/LowQual/giab_lowqual_filtered.bed"\
         "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/benchmarks/metasv/LowQual/metasv_lowqual_filtered.bed")

lowqual_dir="/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standar/HG001/final_gold_standard/general/estudios/lowqual/"
mkdir -p "$lowqual_dir"

# Process each study individually
for input_file in "${inputs[@]}"; do
    filename=$(basename "$input_file")
    cut -f1-4 "$input_file" | bedtools merge -i - -c 4 -o collapse -delim "|" > "$lowqual_dir/$filename"
    echo "Columns (chrom start end svtype) extracted and overlapping variants merged individually for study $filename"
done

# Merge pass and lowqual filtered variants for each study
studies_dir="/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/general/estudios/"
studies=("1000genomes" "giab" "metasv")

for study in "${studies[@]}"; do
    pass_file="$studies_dir/pass/${study}_pass_filtered.bed"
    lowqual_file="$studies_dir/lowqual/${study}_lowqual_filtered.bed"
    # Concatenate sorted lowqual and pass BEDs
    cat "$lowqual_file" "$pass_file" | sort -k1,1V -k2,2n > "$studies_dir/${study}_filtered.bed"
    echo "Pass and lowqual filtered variants are sorted and concatenated for study $study"
done



###################### LAX MERGE FUSION ######################

fusion_dir="/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/general/fusion_merge/"
mkdir -p "$fusion_dir"

# Combine and sort BED files
combined_sorted_bed="${fusion_dir}/estudios_completos_sorted.bed"
cat "${studies_dir}"/*.bed | bedtools sort > "$combined_sorted_bed"

# Merge with variant collapse to include SV type
merged_bed="${fusion_dir}/fusion_merge.bed"
bedtools merge -i "$combined_sorted_bed" -c 4 -o collapse > "$merged_bed"
awk -F'\t' '{if (gsub(",", "&", $4) >= 1) print $0}' OFS='\t' "${merged_bed}" > "${fusion_dir}/fusion_merge_2ormore.bed"
awk -F'\t' '{if (gsub(",", "&", $4) >= 2) print $0}' OFS='\t' "${merged_bed}" > "${fusion_dir}/fusion_merge_3ormore.bed"
awk -F'\t' '{if (gsub(",", "&", $4) >= 3) print $0}' OFS='\t' "${merged_bed}" > "${fusion_dir}/fusion_merge_4.bed"




###################### STRICT INTERSECT FUSION ######################

fusion_dir="/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/general/fusion_intersect/"
mkdir -p "$fusion_dir"

# Multiintersect 
intersect_bed="${fusion_dir}/fusion_intersect.bed"
bedtools  multiinter -i "${studies_dir}"/*.bed -f 0.5 -r > "$intersect_bed"
awk '$4 >= 2' "${fusion_dir}/fusion_intersect.bed" | bedtools merge -i - > "${fusion_dir}/fusion_intersect_2ormore.bed" 
awk '$4 >= 3' "${fusion_dir}/fusion_intersect.bed" | bedtools merge -i - > "${fusion_dir}/fusion_intersect_3ormore.bed"
awk '$4 >= 4' "${fusion_dir}/fusion_intersect.bed" | bedtools merge -i - > "${fusion_dir}/fusion_intersect_4.bed"














