rule exomedepth:
	input:
		bam = os.path.join(config["alignment_dir"], config["sonda"], config["alignment"], "{sample}.bam"),
		ref = config["human_reference"],
		bed = config["exomedepth_bed"],
		controls = lambda wildcards: [os.path.join(config["alignment_dir"], config["sonda"], config["alignment"], control + ".bam") for control in CONTROLS]
		# en el caso de tener un pull sin controles: tomo X muestras del pull y ejecuto X veces tomando 1 muestra de interes y X-1 muestras controles en cada ejecucion
		# controls = lambda wildcards: [os.path.join(config["alignment_dir"], config["sonda"], config["alignment"], f"{ctrl}.bam") for ctrl in SAMPLES if ctrl != wildcards.sample]
	params:
		bam_files = lambda wildcards, input: ','.join([input.bam] + input.controls)
	output:
		tsv =  os.path.join(config["results_dir"], "exomedepth", config["sonda"], config["alignment"], "{sample}", "{sample}.tsv")
	threads: 10
	benchmark: 
		os.path.join(config["results_dir"], "exomedepth", "benchmark", config["sonda"], config["alignment"], "{sample}_exomedepth.log")
	log: 
		os.path.join(config["results_dir"], "exomedepth", "logs", config["sonda"], config["alignment"], "{sample}_exomedepth.log")
	singularity:
		os.path.join(config["singularity_dir"], "exomedepth.sif")
	shell:
		"cd /ingemm/scratch/TFM/CNV/TFM_borrador/algorithms/tools/exomedepth && ./exomedepth_smk.r -b {input.bed} -r {input.ref} -f {params.bam_files} -o {output.tsv}  2> {log}"


rule exomedepth_to_vcf:
	input:
		tsv =  os.path.join(config["results_dir"], "exomedepth", config["sonda"], config["alignment"], "{sample}", "{sample}.tsv"),
	output:
		output_bed = os.path.join(config["results_process_dir"], "exomedepth", config["sonda"], config["alignment"], "{sample}.bed"),
		output_vcf = os.path.join(config["results_process_dir"], "exomedepth", config["sonda"], config["alignment"], "{sample}.vcf"),
	shell:
		"python3 /ingemm/scratch/TFM/CNV/TFM_borrador/algorithms/tools/exomedepth/exomedepth_to_vcf.py --f {input.tsv} --o {output.output_bed}"


rule wittyer_general_exomedepth:
    input: 
        output_vcf = os.path.join(config["results_process_dir"], "exomedepth", config["sonda"], config["alignment"], "{sample}.vcf"),
        gold_standard_general = config["gold_standard_general"]
    output:
        directory(os.path.join(config["results_process_dir"], "exomedepth", config["sonda"], config["alignment"], "{sample}", "general"))
    singularity:
        os.path.join(config["singularity_dir"], "wittyer.sif")
    shell:
        "dotnet /opt/Wittyer/Wittyer.dll -i {input.output_vcf} -t {input.gold_standard_general} -o {output} --pd 1 --bpd 10000 --if PASS --em sc"


rule wittyer_pass_exomedepth:
    input: 
        output_vcf = os.path.join(config["results_process_dir"], "exomedepth", config["sonda"], config["alignment"], "{sample}.vcf"),
        gold_standard_pass =  config["gold_standard_pass"]
    output:
        directory(os.path.join(config["results_process_dir"], "exomedepth", config["sonda"], config["alignment"], "{sample}", "pass"))
    singularity:
        os.path.join(config["singularity_dir"], "wittyer.sif")
    shell:
        "dotnet /opt/Wittyer/Wittyer.dll -i {input.output_vcf} -t {input.gold_standard_pass} -o {output} --pd 1 --bpd 10000 --if PASS --em sc"

		