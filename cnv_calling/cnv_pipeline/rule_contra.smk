rule contra_baseline:
	input:
		bed = config["exome_bed"],
		controls = lambda wildcards: [os.path.join(config["alignment_dir"], config["sonda"], config["alignment"], control + ".bam") for control in CONTROLS],
	params:
		control_bams = lambda wildcards, input: ' '.join(input.controls),
		baseline_path = os.path.join(config["tools_dir"], "contra", "CONTRA.v2.0.8", "baseline.py"),
	output:
		directory(os.path.join(config["results_dir"], "contra", config["sonda"], config["alignment"], "baseline")),
	threads: 10
	singularity:
		os.path.join(config["singularity_dir"], "contra.sif")
	shell:
		"python /contra/baseline.py --target {input.bed} --files {params.control_bams} --output {output} --name baseline_test"
		
rule contra:
	input:
		bam = os.path.join(config["alignment_dir"], config["sonda"], config["alignment"], "{sample}.bam"),
		ref = config["human_reference"],
		bed = config["exome_bed"],
		baseline = os.path.join(config["results_dir"], "contra", config["sonda"], config["alignment"], "baseline", "baseline_test.pooled2_TRIM0.2.txt"),
	output:
		directory(os.path.join(config["results_dir"], "contra", config["sonda"], config["alignment"], "{sample}", "results"))
	threads: 10
	benchmark: 
		os.path.join(config["results_dir"], "contra", "benchmark", config["sonda"], config["alignment"], "{sample}_contra.log")
	log: 
		os.path.join(config["results_dir"], "contra", "logs", config["sonda"], config["alignment"], "{sample}_contra.log")
	singularity:
		os.path.join(config["singularity_dir"], "contra.sif")
	shell:
		"(/contra/contra.py --target {input.bed} --test {input.bam} -c {input.baseline} --bed -o {output}) 2> {log}"


rule contra_to_vcf:
	input:
		txt = os.path.join(config["results_dir"], "contra", config["sonda"], config["alignment"], "{sample}", "results", "table", "CNATable.10rd.10bases.20bins.DetailsFILTERED.txt"),
	output:
		output_bed = os.path.join(config["results_process_dir"], "contra", config["sonda"], config["alignment"], "{sample}.bed"),
		output_vcf = os.path.join(config["results_process_dir"], "contra", config["sonda"], config["alignment"], "{sample}.vcf"),
	shell:
		"python3 /ingemm/scratch/TFM/CNV/TFM_borrador/algorithms/tools/contra/contra_to_vcf.py --f {input.txt} --o {output.output_bed}"


rule wittyer_general_contra:
    input: 
        output_vcf = os.path.join(config["results_process_dir"], "contra", config["sonda"], config["alignment"], "{sample}.vcf"),
        gold_standard_general = config["gold_standard_general"]
    output:
        directory(os.path.join(config["results_process_dir"], "contra", config["sonda"], config["alignment"], "{sample}", "general"))
    singularity:
        os.path.join(config["singularity_dir"], "wittyer.sif")
    shell:
        "dotnet /opt/Wittyer/Wittyer.dll -i {input.output_vcf} -t {input.gold_standard_general} -o {output} --pd 1 --bpd 10000 --if PASS --em sc"


rule wittyer_pass_contra:
    input: 
        output_vcf = os.path.join(config["results_process_dir"], "contra", config["sonda"], config["alignment"], "{sample}.vcf"),
        gold_standard_pass =  config["gold_standard_pass"]
    output:
        directory(os.path.join(config["results_process_dir"], "contra", config["sonda"], config["alignment"], "{sample}", "pass"))
    singularity:
        os.path.join(config["singularity_dir"], "wittyer.sif")
    shell:
        "dotnet /opt/Wittyer/Wittyer.dll -i {input.output_vcf} -t {input.gold_standard_pass} -o {output} --pd 1 --bpd 10000 --if PASS --em sc"

