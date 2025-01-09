import pandas as pd
import os

def df_bed(path):
    cols = ["chr", "start", "end", "svtype"]
    data = pd.read_csv(path, sep="\t", names=cols)
    data["length"] = data["end"] - data["start"]
    sample_name = os.path.basename(path)
    return data, sample_name

dic_gs = {
    "pass" : ["/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/pass/gold_standard_pass_IDTv1.bed",
            "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/pass/gold_standard_pass_IDTv2.bed",
            "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/pass/gold_standard_pass_ROCHEv1.bed"],
            
    "general" : ["/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/general/gold_standard_general_IDTv1.bed",
            "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/general/gold_standard_general_IDTv2.bed",
            "/ingemm/scratch/TFM/CNV/TFM_borrador/gold_standard/HG001/final_gold_standard/general/gold_standard_general_ROCHEv1.bed"]
		}


for path in dic_gs["general"]:
    data, sample_name = df_bed(path)
    print(sample_name)
    print(data["length"].describe())
    print("")
 
 for path in dic_gs["pass"]:
    data, sample_name = df_bed(path)
    print(sample_name)
    print(data["length"].describe())
    print("")
