rule cnvkit:
	input:
		bam = os.path.join(config["alignment_dir"], config["sonda"], config["alignment"], "{sample}.bam"),
		#bai = os.path.join(config["alignment_dir"], config["sonda"], config["alignment"], "{sample}.bam.bai"),
		ref = config["human_reference"],
		bed = config["exome_bed"],
		controls = lambda wildcards: [os.path.join(config["alignment_dir"], config["sonda"], config["alignment"], control + ".bam") for control in CONTROLS],
		refFlat = config["refFlat"],
		mappable = config["bed_mappable"], 
	params:
		outdir = os.path.join(config["results_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}")
	output:
		antitargetcoverage = os.path.join(config["results_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}", "{sample}.antitargetcoverage.cnn"),
		cnr = os.path.join(config["results_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}", "{sample}.cnr"),
		cns = os.path.join(config["results_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}", "{sample}.cns"),
		targetcoverage = os.path.join(config["results_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}", "{sample}.targetcoverage.cnn"),
		diagram = os.path.join(config["results_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}", "{sample}-diagram.pdf"),
		scatter = os.path.join(config["results_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}","{sample}-scatter.pdf"),

		ref_cnn = os.path.join(config["results_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}", "myref.cnn"), 	
	threads: 8
	benchmark: 
		os.path.join(config["results_dir"], "cnvkit", "benchmark", config["sonda"], config["alignment"], "{sample}_cnvkit.log")
	log: 
		os.path.join(config["results_dir"], "cnvkit", "logs", config["sonda"], config["alignment"], "{sample}_cnvkit.log")
	singularity:
		os.path.join(config["singularity_dir"], "cnvkit.sif")
	shell:
		"""
		/usr/bin/cnvkit batch {input.bam} --normal {input.controls} \
		-p {threads} --fasta {input.ref} --annotate {input.refFlat} \
		--access {input.mappable} --method hybrid --targets {input.bed} \
		--output-reference {output.ref_cnn} \
		--output-dir {params.outdir} --scatter --diagram  2> {log}
		"""

rule cnvkit_call:
	input:
		cns = os.path.join(config["results_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}", "{sample}.cns"),
	output:
		call_cns = os.path.join(config["results_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}", "{sample}_call.cns")
	threads: 1
	benchmark: 
		os.path.join(config["results_dir"], "cnvkit", "benchmark", config["sonda"], config["alignment"], "{sample}_cnvkit_call.log")
	log: 
		os.path.join(config["results_dir"], "cnvkit", "logs", config["sonda"], config["alignment"], "{sample}_cnvkit_call.log")
	singularity:
		os.path.join(config["singularity_dir"], "cnvkit.sif")
	shell:
		"/usr/bin/cnvkit call {input.cns} -o {output.call_cns};  2> {log} "


rule cnvkit_bed:
	input:
		call_cns = os.path.join(config["results_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}", "{sample}_call.cns"),
	output:
		call_bed = os.path.join(config["results_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}", "{sample}.bed")
	threads: 1
	benchmark: 
		os.path.join(config["results_dir"], "cnvkit", "benchmark", config["sonda"], config["alignment"], "{sample}_cnvkit_bed.log")
	log: 
		os.path.join(config["results_dir"], "cnvkit", "logs", config["sonda"], config["alignment"], "{sample}_cnvkit_bed.log")
	singularity:
		os.path.join(config["singularity_dir"], "cnvkit.sif")
	shell:
		"/usr/bin/cnvkit export bed {input.call_cns} --show all -y -o {output.call_bed};  2> {log} "


rule cnvkit_to_vcf:
	input:
		call_cns = os.path.join(config["results_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}", "{sample}_call.cns")
	output:
		output_bed = os.path.join(config["results_process_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}.bed"),
		output_vcf = os.path.join(config["results_process_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}.vcf"),
	shell:
		"python3 /ingemm/scratch/TFM/CNV/TFM_borrador/algorithms/tools/cnvkit/cnvkit_to_vcf.py --f {input.call_cns} --o {output.output_bed}"


rule wittyer_general_cnvkit:
    input: 
        output_vcf = os.path.join(config["results_process_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}.vcf"),
        gold_standard_general = config["gold_standard_general"]
    output:
        directory(os.path.join(config["results_process_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}", "general"))
    singularity:
        os.path.join(config["singularity_dir"], "wittyer.sif")
    shell:
        "dotnet /opt/Wittyer/Wittyer.dll -i {input.output_vcf} -t {input.gold_standard_general} -o {output} --pd 1 --bpd 10000 --if PASS --em sc"


rule wittyer_pass_cnvkit:
    input: 
        output_vcf = os.path.join(config["results_process_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}.vcf"),
        gold_standard_pass = config["gold_standard_pass"]
    output:
        directory(os.path.join(config["results_process_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}", "pass"))
    singularity:
        os.path.join(config["singularity_dir"], "wittyer.sif")
    shell:
        "dotnet /opt/Wittyer/Wittyer.dll -i {input.output_vcf} -t {input.gold_standard_pass} -o {output} --pd 1 --bpd 10000 --if PASS --em sc"


