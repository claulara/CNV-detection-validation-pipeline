rule preprocess_bed:
	input:
		bed = config["lacon_exome_bed"],
	output:
		cnv_bed = os.path.join(config["results_dir"], "laconv", config["sonda"], config["alignment"], "preprocess", "240722_Exomas_annotated_for_coverage.preprocessed.bed")
	params:
		cnv_path = os.path.join(config["results_dir"], "laconv", config["sonda"], config["alignment"], "preprocess"),
		pipeline_path = os.path.join(config["tools_dir"], "laconv")
	singularity:
		os.path.join(config["singularity_dir"], "centos_anaconda_pipeline_cnv.sif")
	benchmark: 
		os.path.join(config["results_dir"], "laconv", "benchmark", config["sonda"], config["alignment"], "laconv_preprocess_bed.log")
	log: 
		os.path.join(config["results_dir"], "laconv", "logs", config["sonda"], config["alignment"], "laconv_preprocess_bed.log")
	shell:
		"scif run python_pipeline python {params.pipeline_path}/preprocess_bed_analysis.py --b {input.bed} --o {params.cnv_path}"



rule coverage:
	input:
		bam = os.path.join(config["alignment_dir"], config["sonda"], config["alignment"], "{sample}.bam"),
		cnv_bed = os.path.join(config["results_dir"], "laconv", config["sonda"], config["alignment"], "preprocess", "240722_Exomas_annotated_for_coverage.preprocessed.bed"),
	output:
		regions = os.path.join(config["results_dir"], "laconv", config["sonda"], config["alignment"], "coverage", "{sample}_Mosdepthcoverage.regions.bed"),
	params:
		pipeline_path = os.path.join(config["tools_dir"], "laconv"),
		mosdepth_path = "/scif/apps/mosdepth/bin/mosdepth",
		output_path = os.path.join(config["results_dir"], "laconv", config["sonda"], config["alignment"], "coverage"),
	singularity:
		os.path.join(config["singularity_dir"], "centos_anaconda_pipeline_cnv.sif")
	benchmark: 
		os.path.join(config["results_dir"], "laconv", "benchmark", config["sonda"], config["alignment"], "{sample}_laconv_coverage.log")
	log: 
		os.path.join(config["results_dir"], "laconv", "logs", config["sonda"], config["alignment"], "{sample}_laconv_coverage.log")
	shell:
		"scif run python_pipeline python {params.pipeline_path}/run_coverage.py --b {input.cnv_bed} --m {params.mosdepth_path} --o {params.output_path} --f {input.bam}"
		


rule lacon_cnvs:
	input:
		cfg = config["lacon_config"],
		bed = os.path.join(config["results_dir"], "laconv", config["sonda"], config["alignment"], "preprocess", "240722_Exomas_annotated_for_coverage.preprocessed.bed"),
		regions = expand(os.path.join(config["results_dir"], "laconv", config["sonda"], config["alignment"], "coverage", "{sample}_Mosdepthcoverage.regions.bed"), sample=ALL_SAMPLES), 
		bam = expand(os.path.join(config["alignment_dir"], config["sonda"], config["alignment"], "{sample}.bam"), sample=ALL_SAMPLES), 
		bai = expand(os.path.join(config["alignment_dir"], config["sonda"], config["alignment"], "{sample}.bam.bai"), sample=ALL_SAMPLES), 
	output:
		laconv_bed =  expand(os.path.join(config["results_dir"], "laconv", config["sonda"], config["alignment"], "laconv", "{sample}.cnv.bed"), sample=ALL_SAMPLES),		   
	params:
		cnv_path = os.path.join(config["results_dir"], "laconv", config["sonda"], config["alignment"]),
		pipeline_path =  os.path.join(config["tools_dir"], "laconv"),
		p_bam = ("--f " + str.join(" --f ", expand(os.path.join(config["alignment_dir"], config["sonda"], config["alignment"], "{sample}.bam"), zip, sample=ALL_SAMPLES))),
	singularity:
		os.path.join(config["singularity_dir"], "centos_anaconda_pipeline_cnv.sif"),
	benchmark: 
		os.path.join(config["results_dir"], "laconv", "benchmark", config["sonda"], config["alignment"], "laconv_run.log")
	log: 
		os.path.join(config["results_dir"], "laconv", "logs", config["sonda"], config["alignment"], "laconv_run.log")
	shell:
		"scif run python_pipeline python {params.pipeline_path}/CNV_analysis.py --no-cov --clean --cfg {input.cfg} --b {input.bed} {params.p_bam} --o {params.cnv_path}"



rule laconv_to_vcf:
	input:
		laconv_bed = os.path.join(config["results_dir"], "laconv", config["sonda"], config["alignment"], "{sample}", "laconv", "{sample}.cnv.bed"),
	output:
		output_bed = os.path.join(config["results_process_dir"], "laconv", config["sonda"], config["alignment"], "{sample}.bed"),
		output_vcf = os.path.join(config["results_process_dir"], "laconv", config["sonda"], config["alignment"], "{sample}.vcf"),
	shell:
		"python3 /ingemm/scratch/TFM/CNV/TFM_borrador/algorithms/tools/laconv/laconv_to_vcf.py --f {input.laconv_bed} --o {output.output_bed}"


rule wittyer_general_laconv:
    input: 
        output_vcf = os.path.join(config["results_process_dir"], "laconv", config["sonda"], config["alignment"], "{sample}.vcf"),
        gold_standard_general = config["gold_standard_general"]
    output:
        directory(os.path.join(config["results_process_dir"], "laconv", config["sonda"], config["alignment"], "{sample}", "general"))
    singularity:
        os.path.join(config["singularity_dir"], "wittyer.sif")
    shell:
        "dotnet /opt/Wittyer/Wittyer.dll -i {input.output_vcf} -t {input.gold_standard_general} -o {output} --pd 1 --bpd 10000 --if PASS --em sc"


rule wittyer_pass_laconv:
    input: 
        output_vcf = os.path.join(config["results_process_dir"], "laconv", config["sonda"], config["alignment"], "{sample}.vcf"),
        gold_standard_pass =  config["gold_standard_pass"]
    output:
        directory(os.path.join(config["results_process_dir"], "laconv", config["sonda"], config["alignment"], "{sample}", "pass"))
    singularity:
        os.path.join(config["singularity_dir"], "wittyer.sif")
    shell:
        "dotnet /opt/Wittyer/Wittyer.dll -i {input.output_vcf} -t {input.gold_standard_pass} -o {output} --pd 1 --bpd 10000 --if PASS --em sc"


