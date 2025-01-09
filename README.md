
# Validation set and Snakemake Workflow for CNV detection on NGS data

Copy Number Variations (CNV) are one of the most significant contributors to genetic variation in the human genome and play a crucial role in various diseases and phenotypes. However, their detection and characterization remain challenging due to intrinsic limitations of sequencing technologies and genomic complexities, such as repetitive or complex regions that are difficult to capture.

Therefore, this study proposes the following objectives:

 - Develop a robust and reliable validation set of CNVs for the NA12878 sample, a widely studied reference in the scientific community.
 - Implement 7 bioinformatics algorithms for CNV detection within an optimized workflow, using Singularity and Snakemake.
 - Evaluate the performance of these algorithms against the previously generated validation set with precision, recall, and f-score metrics.

## Validation Set for the NA12878 Sample

To generate a robust validation set or gold standard, we used the following four public repositories, that are widely recognized and used within the scientific community. These repositories provide structural variant datasets specifically for the NA12878 sample.

```
# The 1000 Genome Project 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ 
# Genome in a Bottle (GIAB)
https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/
# svclassify (Parikh et al.)
https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/
https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/metasv_trio_validation/

```
The development of this reference set is located in the directory *gold_standard*, organized with the following folder structure:

```
gold_standard
├── gold_standard
│   ├── goldstandard_statistics.py           # Python script to generate statistics of the gold standard
│   ├── goldstandard_to_vcf.py               # Python script to convert gold standard to VCF format for witty.er
│   ├── general                              # General gold standard
│   │   ├── create_goldstandard_general.sh   # Shell script to create the general gold standard
│   │   ├── fusion_intersect                 # Bedtools multiinter fusion results
│   │   ├── fusion_merge                     # Bedtools merge fusion results
│   └── pass                                 # Pass gold standard
│       ├── create_goldstandard_pass.sh      # Shell script to create the pass gold standard
│       ├── fusion_intersect                 # Bedtools multiinter fusion results
│       ├── fusion_merge                     # Bedtools merge fusion results
└── studies_resources                        # Resources and data for the studies
    ├── 1000genomes                          # 1000 Genomes data
    │   ├── LowQual                          # Low-quality variants
    │   └── PASS                             # High-confidence variants (PASS)
    ├── giab                                 # Genome in a Bottle (GIAB) data
    │   ├── LowQual                          # Low-quality variants
    │   └── PASS                             # High-confidence variants (PASS)
    ├── metasv                               # MetaSV data
    │   ├── LowQual                          # Low-quality variants
    │   └── PASS                             # High-confidence variants (PASS)
    ├── svclassify                           # SVClassify data
    ├── studies_process_to_bed.py            # Python script to process study data into BED format
    └── studies_statistics_plots.py          # Python script to generate plots and statistics of each study
```


## Implementation and Validation of a Bioinformatics Pipeline for CNV Detection
The CNV detection algorithms implemented in this study are listed below along with their corresponding repositories:
- [cn.MOPS](https://bioconductor.org/packages/release/bioc/html/cn.mops.html)
- [CNVkit](https://cnvkit.readthedocs.io/en/stable/)
- [CONTRA](https://contra-cnv.sourceforge.net/)
- [ExomeDepth](https://github.com/vplagnol/ExomeDepth)
- [Manta](https://github.com/Illumina/manta)
- LACONv (an internally developed algorithm by the Genetics Service at INGEMM, La Paz University Hospital)
- [XHMM](https://statgen.bitbucket.io/xhmm/index.html)



To develop a unified workflow that integrates all these algorithms into a single process, ensuring reproducibility and optimization for each, the ```Snakemake``` workflow engine is employed.  his tool automates and scales complex workflows through a structured set of rules, where each rule represents a step in the process and defines the inputs, outputs, and commands required for execution.

Therefore, each of the aforementioned algorithms is implemented in the  ```Snakemake``` language, creating a separate *Snakefile* for each algorithm. These *Snakefiles* define the required rules to execute the steps of each algorithm.

Finally, a general Snakefile named snakefile.smk is created, which includes all these rules and runs the entire pipeline from start to finish with a single command. The command to execute this pipeline is as follows:

```
snakemake -j 20 -s /ingemm/scratch/TFM/CNV/TFM_borrador/algorithms/pipeline/snakefile.smk 
          --profile /ingemm/scratch/TFM/CNV/TFM_borrador/algorithms/profile/ 
          --configfile /ingemm/scratch/TFM/CNV/TFM_borrador/algorithms/pipeline/configuraciones/IDT-V1/config_b37_minimap.yaml
```

The data in this study were sequenced using different exome capture kits and processed through various bioinformatics workflows, combining different aligners and read filtering tools. 9 configuration files were created to represent all possible scenarios. So the pipeline was executed using the previous command for each configuration.

The configuration files are located in the *configurations_yaml* directory, along with the corresponding groups of samples and controls associated with each configuration.

The validation was performed using the ```witty.er``` algorithm, which is a tool designed to evaluate the precision, recall, and f-score of variant detection methods by comparing predicted variants to a reference set.

- [witty.er](https://github.com/Illumina/witty.er)


The code required to implement each CNV detection algorithm, their integration into the ```Snakemake``` workflow, and the results obtained from the validation against the generated gold standards using ```witty.er```, are organized with the following folder structure:
```
cnv_calling                                           # Main directory for CNV calling analysis
├── algorithm_results                                 # Results and validation of CNV calling algorithms
│   ├── cnv_calling_results_normalized                # Normalized CNV calling results for each algorithm on BED and VCF format
│   ├── plotting_scripts                              # Python scripts to generate plots for algorithm evaluation
│   │   ├── algorithm_performance_plots.py            # Plots comparing overall algorithm performance
│   │   ├── heatmap_global_comparation.py             # Heatmaps for global comparisons between algorithms
│   │   ├── kit_sample_pipeline_performance_plots.py  # Visualizations of kit, sample, and pipeline performance
│   │   └── svtype_performance_plots.py               # Plots for performance analysis by structural variant type
│   └── witty_results                                 # Witty.er tool results and metrics
│       ├── json_to_df.py                             # Script to parse JSON witty.er results into dataframes
│       ├── overall_metrics.csv                       # CSV file with overall metrics from algorithms
│       └── svtype_metrics.csv                        # CSV file with metrics by structural variant type 
├── cnv_algorithms                                    # CNV detection algorithms
│   ├── cnmops                                        # CNMOPS algorithm resources
│   ├── cnvkit                                        # CNVkit algorithm resources
│   ├── contra                                        # CONTRA algorithm resources
│   ├── exomedepth                                    # ExomeDepth algorithm resources
│   ├── laconv                                        # LACONv algorithm resources
│   ├── manta                                         # Manta algorithm resources
│   └── xhmm                                          # XHMM algorithm directory
└── cnv_pipeline                                      # Snakemake workflow for CNV calling: rules and configurations
   └── configurations_yaml                            # Configuration files for pipeline and capture kit settings

```
