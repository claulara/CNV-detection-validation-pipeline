rule manta:
	input:
		bam = os.path.join(config["alignment_dir"], config["sonda"], config["alignment"], "{sample}.bam"),
		ref = config["human_reference"]
	output:
		vcf_diploid = os.path.join(config["results_dir"], "manta", config["sonda"], config["alignment"], "{sample}", "results", "variants", "diploidSV.vcf.gz"),
		tbi_diploid = os.path.join(config["results_dir"], "manta", config["sonda"], config["alignment"], "{sample}", "results", "variants", "diploidSV.vcf.gz.tbi"),
		vcf_small = os.path.join(config["results_dir"], "manta", config["sonda"], config["alignment"], "{sample}", "results", "variants", "candidateSmallIndels.vcf.gz"),
		tbi_small = os.path.join(config["results_dir"], "manta", config["sonda"], config["alignment"], "{sample}", "results", "variants", "candidateSmallIndels.vcf.gz.tbi"),
		vcf_cand = os.path.join(config["results_dir"], "manta", config["sonda"], config["alignment"], "{sample}", "results", "variants", "candidateSV.vcf.gz"),
		tbi_cand = os.path.join(config["results_dir"], "manta", config["sonda"], config["alignment"], "{sample}", "results", "variants", "candidateSV.vcf.gz.tbi"),
		stats_summary = os.path.join(config["results_dir"], "manta", config["sonda"], config["alignment"], "{sample}", "results", "stats", "alignmentStatsSummary.txt"),
		stats_cand_tsv = os.path.join(config["results_dir"], "manta", config["sonda"], config["alignment"], "{sample}", "results", "stats", "svCandidateGenerationStats.tsv"),
		stats_cand_xml = os.path.join(config["results_dir"], "manta", config["sonda"], config["alignment"], "{sample}", "results", "stats", "svCandidateGenerationStats.xml"),
		stats_graph = os.path.join(config["results_dir"], "manta", config["sonda"], config["alignment"], "{sample}", "results", "stats", "svLocusGraphStats.tsv"),
	params:
		outdir = os.path.join(config["results_dir"], "manta", config["sonda"], config["alignment"], "{sample}")
	threads: 1
	benchmark: 
		os.path.join(config["results_dir"], "manta", "benchmark", config["sonda"], config["alignment"], "{sample}_manta.log")
	log: 
		os.path.join(config["results_dir"], "manta", "logs", config["sonda"], config["alignment"], "{sample}_manta.log")
	singularity:
		os.path.join(config["singularity_dir"], "manta.sif")
	shell:
		# Configure Manta
		"/usr/bin/manta/bin/configManta.py  "
		"--bam {input.bam} "
		"--referenceFasta {input.ref} "
		"--runDir {params.outdir} "
		"--exome;"
		# Run Manta
		"{params.outdir}/runWorkflow.py --mode local --jobs {threads}  2> {log}"
		

rule manta_to_vcf:
	input:
		vcf_diploid = os.path.join(config["results_dir"], "manta", config["sonda"], config["alignment"], "{sample}", "results", "variants", "diploidSV.vcf.gz"),
	output:
		output_bed = os.path.join(config["results_process_dir"], "manta", config["sonda"], config["alignment"], "{sample}.bed"),
		output_vcf = os.path.join(config["results_process_dir"], "manta", config["sonda"], config["alignment"], "{sample}.vcf"),
	shell:
		"python3 /ingemm/scratch/TFM/CNV/TFM_borrador/algorithms/tools/manta/manta_to_vcf.py --f {input.vcf_diploid} --o {output.output_bed}"


rule wittyer_general_manta:
    input: 
        output_vcf = os.path.join(config["results_process_dir"], "manta", config["sonda"], config["alignment"], "{sample}.vcf"),
        gold_standard_general = config["gold_standard_general"]
    output:
        directory(os.path.join(config["results_process_dir"], "manta", config["sonda"], config["alignment"], "{sample}", "general"))
    singularity:
        os.path.join(config["singularity_dir"], "wittyer.sif")
    shell:
        "dotnet /opt/Wittyer/Wittyer.dll -i {input.output_vcf} -t {input.gold_standard_general} -o {output} --pd 1 --bpd 10000 --if PASS --em sc"


rule wittyer_pass_manta:
    input: 
        output_vcf = os.path.join(config["results_process_dir"], "manta", config["sonda"], config["alignment"], "{sample}.vcf"),
        gold_standard_pass =  config["gold_standard_pass"]
    output:
        directory(os.path.join(config["results_process_dir"], "manta", config["sonda"], config["alignment"], "{sample}", "pass"))
    singularity:
        os.path.join(config["singularity_dir"], "wittyer.sif")
    shell:
        "dotnet /opt/Wittyer/Wittyer.dll -i {input.output_vcf} -t {input.gold_standard_pass} -o {output} --pd 1 --bpd 10000 --if PASS --em sc"


