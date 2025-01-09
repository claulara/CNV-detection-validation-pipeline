shell.prefix("umask g+w; ")

import os
import pandas as pd

# Load samples
samples_table=pd.read_table(config['ss_samples'], index_col=None)
SAMPLES=samples_table['Samples'].tolist()

# Load controls
controls_table=pd.read_table(config['ss_controls'], index_col=None)
CONTROLS=controls_table['Samples'].tolist()

pipeline_dir = config['pipeline_dir']
include: os.path.join(pipeline_dir, 'rule_cnmops.smk')
include: os.path.join(pipeline_dir, 'rule_cnvkit.smk')
include: os.path.join(pipeline_dir, 'rule_contra.smk')
include: os.path.join(pipeline_dir, 'rule_manta.smk')
include: os.path.join(pipeline_dir, 'rule_lacon.smk')
include: os.path.join(pipeline_dir, 'rule_xhmm.smk')

rule all: 
	input:
		# cnmops
		expand(os.path.join(config["results_process_dir"], "cnmops", config["sonda"], config["alignment"], "{sample}", "general"), sample=SAMPLES),
		expand(os.path.join(config["results_process_dir"], "cnmops", config["sonda"], config["alignment"], "{sample}", "pass"), sample=SAMPLES),
		# cnvkit
		expand(os.path.join(config["results_process_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}", "general"), sample=SAMPLES),
		expand(os.path.join(config["results_process_dir"], "cnvkit", config["sonda"], config["alignment"], "{sample}", "pass"), sample=SAMPLES),
		# contra
		expand(os.path.join(config["results_process_dir"], "contra", config["sonda"], config["alignment"], "{sample}", "general"), sample=SAMPLES),
		expand(os.path.join(config["results_process_dir"], "contra", config["sonda"], config["alignment"], "{sample}", "pass"), sample=SAMPLES),
		# exomedepth
		expand(os.path.join(config["results_process_dir"], "exomedepth", config["sonda"], config["alignment"], "{sample}", "general"), sample=SAMPLES),
		expand(os.path.join(config["results_process_dir"], "exomedepth", config["sonda"], config["alignment"], "{sample}", "pass"), sample=SAMPLES),
		# manta
		expand(os.path.join(config["results_process_dir"], "manta", config["sonda"], config["alignment"], "{sample}", "general"), sample=SAMPLES),
		expand(os.path.join(config["results_process_dir"], "manta", config["sonda"], config["alignment"], "{sample}", "pass"), sample=SAMPLES),
		# laconv
		expand(os.path.join(config["results_process_dir"], "laconv", config["sonda"], config["alignment"], "{sample}", "general"), sample=SAMPLES),
		expand(os.path.join(config["results_process_dir"], "laconv", config["sonda"], config["alignment"], "{sample}", "pass"), sample=SAMPLES),
		# xhmm
		expand(os.path.join(config["results_process_dir"], "xhmm", config["sonda"], config["alignment"], "{sample}", "general"), sample=SAMPLES),
		expand(os.path.join(config["results_process_dir"], "xhmm", config["sonda"], config["alignment"], "{sample}", "pass"), sample=SAMPLES),
		