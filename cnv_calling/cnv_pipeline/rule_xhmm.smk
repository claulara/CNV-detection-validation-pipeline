rule mosdepth:
	input:
		bam = os.path.join(config["alignment_dir"], config["sonda"], config["alignment"], "{sample}.bam"),
		bai = os.path.join(config["alignment_dir"], config["sonda"], config["alignment"], "{sample}.bam.bai"),
		bed = config["exome_bed"],
	params:
		MQ = "15",
		sample_name = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "{sample}", "mosdepth", "{sample}"),
	output:
		temp(os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "{sample}", "mosdepth", "{sample}.mosdepth.global.dist.txt")),
		temp(os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "{sample}", "mosdepth", "{sample}.mosdepth.region.dist.txt")),
		temp(os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "{sample}", "mosdepth", "{sample}.mosdepth.summary.txt")),
		temp(os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "{sample}", "mosdepth", "{sample}.regions.bed.gz")),
		temp(os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "{sample}", "mosdepth", "{sample}.regions.bed.gz.csi")),
	threads: 10
	benchmark: 
		os.path.join(config["results_dir"], "xhmm", "benchmark", config["sonda"], config["alignment"], "{sample}_xhmm_mosdepth.log")
	log: 
		os.path.join(config["results_dir"], "xhmm", "logs", config["sonda"], config["alignment"], "{sample}_xhmm_mosdepth.log")
	singularity:
		os.path.join(config["singularity_dir"], "mosdepth.sif")
	shell:
		"/usr/bin/mosdepth -n -Q {params.MQ} -t {threads} --by {input.bed} {params.sample_name} {input.bam} 2> {log}"


rule merge_coverage:
	input:
		expand(os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "{sample}", "mosdepth", "{sample}.regions.bed.gz"), sample = ALL_SAMPLES),
	params:
		all_regions_bed = lambda wildcards, input: ','.join([f for f in input]),
		mat_cov_py_path = os.path.join(config["tools_dir"], "xhmm", "merge_coverage.py"),
	benchmark: 
		os.path.join(config["results_dir"], "xhmm", "benchmark", config["sonda"], config["alignment"], "xhmm_merge_coverage.log")
	log: 
		os.path.join(config["results_dir"], "xhmm", "logs", config["sonda"], config["alignment"], "xhmm_merge_coverage.log")
	output: 
		matrix_cov = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "coverage_matrix.tsv"),
	shell:
		"python3 {params.mat_cov_py_path} --b {params.all_regions_bed} --o {output} 2> {log}"


rule center_filter:
	input:
		matrix_cov = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "coverage_matrix.tsv"),
	params:
		xhmm_path = os.path.join("/scif/apps/xhmm/bin/xhmm"),
		minTargetSize = "10",
		maxTargetSize = "10000",
		minMeanTargetRD = "10",
		maxMeanTargetRD = "2000",
		minMeanSampleRD = "10",
		maxMeanSampleRD = "2000",
		minSdSampleRD = "0",
		maxSdSampleRD = "500",
	output:
		matrix_cov_filtered_center = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "coverage_matrix_filtered_center.tsv"),
		matrix_cov_excluded_targets = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "coverage_matrix_filtered_center_excluded_targets.tsv"),
		matrix_cov_excluded_samples = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "coverage_matrix_filtered_center_excluded_samples.tsv"),
	benchmark: 
		os.path.join(config["results_dir"], "xhmm", "benchmark", config["sonda"], config["alignment"], "xhmm_center_filter.log")
	log: 
		os.path.join(config["results_dir"], "xhmm", "logs", config["sonda"], config["alignment"], "xhmm_center_filter.log")
	singularity:
		os.path.join(config["singularity_dir"], "centos_anaconda_pipeline_cnv.sif"),
	shell:
		" {params.xhmm_path} --matrix --centerData --centerType target -r {input.matrix_cov} -o {output.matrix_cov_filtered_center}  \
		--outputExcludedTargets {output.matrix_cov_excluded_targets} \
		--outputExcludedSamples {output.matrix_cov_excluded_samples}  \
		--minTargetSize {params.minTargetSize} --maxTargetSize {params.maxTargetSize} \
		--minMeanTargetRD {params.minMeanTargetRD} --maxMeanTargetRD {params.maxMeanTargetRD} \
		--minMeanSampleRD {params.minMeanSampleRD} --maxMeanSampleRD {params.maxMeanSampleRD} \
		--minSdSampleRD {params.minSdSampleRD} --maxSdSampleRD {params.maxSdSampleRD} 2> {log}"


rule pca:
	input:
		matrix_cov_filtered_center = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "coverage_matrix_filtered_center.tsv"),
	params:
		xhmm_path = os.path.join("/scif/apps/xhmm/bin/xhmm"),
	output:
		pc = directory(os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "pca")),
	singularity:
		os.path.join(config["singularity_dir"], "centos_anaconda_pipeline_cnv.sif"),
	benchmark: 
		os.path.join(config["results_dir"], "xhmm", "benchmark", config["sonda"], config["alignment"], "xhmm_pca.log")
	log: 
		os.path.join(config["results_dir"], "xhmm", "logs", config["sonda"], config["alignment"], "xhmm_pca.log")
	shell:
		"""
		mkdir -p {output.pc}
		{params.xhmm_path} --PCA -r {input.matrix_cov_filtered_center} --PCAfiles {output.pc} 2> {log}
		"""


rule normalize:
	input:
		matrix_cov_filtered_center = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "coverage_matrix_filtered_center.tsv"),
		pc = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "pca"),
	params:
		xhmm_path = os.path.join("/scif/apps/xhmm/bin/xhmm"),
	output:
		pca_normalized = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "PCA_normalized.txt"),
	singularity:
		os.path.join(config["singularity_dir"], "centos_anaconda_pipeline_cnv.sif"),
	benchmark: 
		os.path.join(config["results_dir"], "xhmm", "benchmark", config["sonda"], config["alignment"], "xhmm_normalize.log")
	log: 
		os.path.join(config["results_dir"], "xhmm", "logs", config["sonda"], config["alignment"], "xhmm_normalize.log")
	shell:
		" {params.xhmm_path} --normalize -r {input.matrix_cov_filtered_center} \
		 --PCAfiles {input.pc} --normalizeOutput {output.pca_normalized} --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7 2> {log}"


rule normalize_center_filter:
	input:
		pca_normalized = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "PCA_normalized.txt"),
	params:
		xhmm_path = os.path.join("/scif/apps/xhmm/bin/xhmm"),
		maxSdTargetRD = 30
	output:
		pca_norm_filtered_center = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "pca_normalized_zscore_filtered_center.tsv"),
		pca_norm_excluded_targets = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "pca_normalized_zscore_filtered_center_excluded_targets.tsv"),
		pca_norm_excluded_samples = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "pca_normalized_zscore_filtered_center_excluded_samples.tsv"),
	singularity:
		os.path.join(config["singularity_dir"], "centos_anaconda_pipeline_cnv.sif"),
	benchmark: 
		os.path.join(config["results_dir"], "xhmm", "benchmark", config["sonda"], config["alignment"], "xhmm_norm_center_filter.log")
	log: 
		os.path.join(config["results_dir"], "xhmm", "logs", config["sonda"], config["alignment"], "xhmm_norm_center_filter.log")
	shell:
		" {params.xhmm_path} --matrix --centerData --centerType sample --zScoreData -r {input.pca_normalized} -o {output.pca_norm_filtered_center}  \
		--outputExcludedTargets {output.pca_norm_excluded_targets} \
		--outputExcludedSamples {output.pca_norm_excluded_samples}  \
		--maxSdSampleRD 30 2> {log} "


rule filter_original_matrix_cov:
	input:
		matrix_cov = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "coverage_matrix.tsv"),
		matrix_cov_excluded_targets = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "coverage_matrix_filtered_center_excluded_targets.tsv"),
		matrix_cov_excluded_samples = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "coverage_matrix_filtered_center_excluded_samples.tsv"),
		pca_norm_excluded_targets = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "pca_normalized_zscore_filtered_center_excluded_targets.tsv"),
		pca_norm_excluded_samples = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "pca_normalized_zscore_filtered_center_excluded_samples.tsv"),
	params:
		xhmm_path = os.path.join("/scif/apps/xhmm/bin/xhmm"),
	output:
		matrix_cov_filtered = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "final_coverage_matrix_filtered.tsv"),
	singularity:
		os.path.join(config["singularity_dir"], "centos_anaconda_pipeline_cnv.sif"),
	benchmark: 
		os.path.join(config["results_dir"], "xhmm", "benchmark", config["sonda"], config["alignment"], "xhmm_filter_original_mat.log")
	log: 
		os.path.join(config["results_dir"], "xhmm", "logs", config["sonda"], config["alignment"], "xhmm_filter_original_mat.log")
	shell:
		" {params.xhmm_path} --matrix -r {input.matrix_cov} \
		--excludeTargets {input.matrix_cov_excluded_targets} --excludeTargets {input.pca_norm_excluded_targets} \
		--excludeSamples {input.matrix_cov_excluded_samples} --excludeSamples {input.pca_norm_excluded_samples}  \
		-o {output.matrix_cov_filtered} 2> {log}" 


rule xhmm:
	input:
		pca_norm_filtered_center = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "pca_normalized_zscore_filtered_center.tsv"),
		matrix_cov_filtered = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "final_coverage_matrix_filtered.tsv"),
	params:
		xhmm_path = os.path.join("/scif/apps/xhmm/bin/xhmm"),
		posterior_path = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "posterior"),
		output_path = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"]),
	output:
		cnv_output = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "cnv_output.txt"),
		cnv_output_aux = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "cnv_output_aux.txt"),
	singularity:
		os.path.join(config["singularity_dir"], "centos_anaconda_pipeline_cnv.sif"),
	benchmark: 
		os.path.join(config["results_dir"], "xhmm", "benchmark", config["sonda"], config["alignment"], "xhmm_run.log")
	log: 
		os.path.join(config["results_dir"], "xhmm", "logs", config["sonda"], config["alignment"], "xhmm_run.log")
	shell:
		"""
		# crear el archivo de params
		echo "1e-8   6   70  -3  1.00    0   1.00    3   1.00" > {params.output_path}/params.txt
		
		{params.xhmm_path} --discover -p {params.output_path}/params.txt \
		-r {input.pca_norm_filtered_center} -R {input.matrix_cov_filtered} \
		-c {output.cnv_output} -a {output.cnv_output_aux} -s {params.posterior_path} 2> {log}
		"""


rule split_cnvs_sample:
	input:
		cnv_output=os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "cnv_output.txt"),
	params:
		split_sample_py_path=os.path.join(config["tools_dir"], "xhmm", "split_results_cnvs.py"),
	output:
		directory(os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "results")),
	shell:
		"python3 {params.split_sample_py_path} --f {input.cnv_output} --o {output}"

rule xhmm_to_vcf:
	input:
		csv = os.path.join(config["results_dir"], "xhmm", config["sonda"], config["alignment"], "results", "{sample}", "{sample}_cnv_output.csv"),
	output:
		output_bed = os.path.join(config["results_process_dir"], "xhmm", config["sonda"], config["alignment"], "{sample}.bed"),
		output_vcf = os.path.join(config["results_process_dir"], "xhmm", config["sonda"], config["alignment"], "{sample}.vcf"),
	shell:
		"python3 /ingemm/scratch/TFM/CNV/TFM_borrador/algorithms/tools/xhmm/xhmm_to_vcf.py --f {input.csv} --o {output.output_bed}"


rule wittyer_general_xhmm:
    input: 
        output_vcf = os.path.join(config["results_process_dir"], "xhmm", config["sonda"], config["alignment"], "{sample}.vcf"),
        gold_standard_general = config["gold_standard_general"]
    output:
        directory(os.path.join(config["results_process_dir"], "xhmm", config["sonda"], config["alignment"], "{sample}", "general"))
    singularity:
        os.path.join(config["singularity_dir"], "wittyer.sif")
    shell:
        "dotnet /opt/Wittyer/Wittyer.dll -i {input.output_vcf} -t {input.gold_standard_general} -o {output} --pd 1 --bpd 10000 --if PASS --em sc"


rule wittyer_pass_xhmm:
    input: 
        output_vcf = os.path.join(config["results_process_dir"], "xhmm", config["sonda"], config["alignment"], "{sample}.vcf"),
        gold_standard_pass =  config["gold_standard_pass"]
    output:
        directory(os.path.join(config["results_process_dir"], "xhmm", config["sonda"], config["alignment"], "{sample}", "pass"))
    singularity:
        os.path.join(config["singularity_dir"], "wittyer.sif")
    shell:
        "dotnet /opt/Wittyer/Wittyer.dll -i {input.output_vcf} -t {input.gold_standard_pass} -o {output} --pd 1 --bpd 10000 --if PASS --em sc"





