[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17302040.svg)](https://doi.org/10.5281/zenodo.17302040)
# Element vs Illumina benchmark

Project: SR_23_02

## Data organisation

- **`data/wgs`**: Raw FASTQs organized under here. To download data, see [Data download](#data-download) section.
- **`analysis`**: Analysis runs on data, either using Snakemake script or existing Nextflow/nf-core workflow
- **`resources`**: Genome annotations and reference
- **`notebooks`**: Jupyter notebooks for data processing and visualization
- **`scripts`**: Custom Python scripts
- **`figures/svg`**: SVG figures
- **`env`**: Environment files 

### Code to figure/table

Notebooks used to generate each figure/table is specified below:

**Figure 1**
- b: [duplicates.ipynb](notebooks/duplicates.ipynb)
- c: [samtools_stats_all.ipynb](notebooks/samtools_stats_all.ipynb)
- d: [samtools_stats_all.ipynb](notebooks/samtools_stats_all.ipynb)
- e: [differential_coverage.ipynb](notebooks/differential_coverage.ipynb)
- f: [variant_calling_benchmarks.ipynb](notebooks/variant_calling_benchmarks.ipynb)
- g: [variant_calling_benchmarks.ipynb](notebooks/variant_calling_benchmarks.ipynb)
- h: [variant_calling_benchmarks.ipynb](notebooks/variant_calling_benchmarks.ipynb)

**Figure 2**
- a: [samtools_stats_per_read.ipynb](notebooks/samtools_stats_per_read.ipynb) 
- b: [fragurancy_error_rate.ipynb](notebooks/fragurancy_error_rate.ipynb)
- c: [samtools_stats_per_read_insert_size.ipynb](notebooks/samtools_stats_per_read_insert_size.ipynb)
- d: [compare_read_stack_multiple.ipynb](notebooks/compare_read_stack_multiple.ipynb)

**Figure 3**
- b: [g4_soft_clipped.ipynb](notebooks/g4_soft_clipped.ipynb)
- c: [stratification_error_rate.ipynb](notebooks/stratification_error_rate.ipynb)
- d: [compare_read_stack_multiple.ipynb](notebooks/compare_read_stack_multiple.ipynb)
- e: [stratification_error_rate.ipynb](notebooks/stratification_error_rate.ipynb)

**Supplementary Table 5**
[g4_overlap.ipynb](notebooks/g4_overlap.ipynb)

**Supplementary Figure 2**
[samtools_stats_all.ipynb](notebooks/samtools_stats_all.ipynb)

**Supplementary Figure 3**
[duplicates.ipynb](notebooks/duplicates.ipynb)

**Supplementary Figure 4**
[samtools_stats_all.ipynb](notebooks/samtools_stats_all.ipynb)

**Supplementary Figure 5**
[samtools_stats_per_read.ipynb](notebooks/samtools_stats_per_read.ipynb) 

**Supplementary Figure 6**
[samtools_stats_per_read_public.ipynb](notebooks/samtools_stats_per_read_public.ipynb)

**Supplementary Figure 7**
[fragurancy_error_rate.ipynb](notebooks/fragurancy_error_rate.ipynb)

## Reproducing analysis

To reproduce the analysis done for this work perform the following steps

0. Clone (`git clone ...`) this repository 
1. Download [Element/Illumina FASTQs](#illumina--element-fastq-download) and [PacBio BAMs](#pacbio-bam-download)
2. Prepare enviroments and software to run analysis, see [Environment setup and configuration](#environment-setup-and-configuration)
3. Download and prepare [resource files](#resources) e.g. genome. 
4. Run [nf-core/sarek](#running-nf-coresarek-on-elementillumina-data)
5. Run [secondary analysis using snakemake worflows](#running-secondary-snakemake-analysis) 

## Data download

### Illumina + Element FASTQ download

FASTQs are publically available on ENA: https://www.ebi.ac.uk/ena/browser/view/PRJEB90663

To download the short read fastqs there is a bash script `download_fastqs.sh`. I.e run:

```
bash download_fastqs.sh
```

The script will download FASTQs in the same folder structure to be able to run nf-core sarek.

### PacBio BAM download

Aligned BAMs are available on ENA: https://www.ebi.ac.uk/ena/browser/view/PRJEB95775 

There is a download bash script in `data/wgs/PacBio_HiFi_BAMs/download_pacbio_bam_ena.sh`

Run this within the directory i.e.

```
cd data/wgs/PacBio_HiFi_BAMs
bash download_pacbio_bam_ena.sh
```

# Public Element FreeStyle dataset

FASTQs be downloaded from here: https://go.elementbiosciences.com/human-whole-genome-sequencing-third-party (dataset `JM-L825-HG002`). Place the FASTQs in the `data/wgs/Element_Freestyle`
 folder.

## Environment setup and configuration

### Workflow dependencies

The **Snakemake** workflows require `snakemake` v8.20.1 (later version may also work) to run, see [documentation here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

The **nf-core/sarek pipeline** is based on **Nextflow**, see [documentation to setup Nextflow](https://nf-co.re/docs/usage/installation). Nextflow v24.10.4 was used in this work.


### Executables

Some of the Snakemake workflows requires executables to be downloaded and added to the configurations. 

Download the required executables below

- mosdepth_d4: https://github.com/brentp/mosdepth/releases/download/v0.3.10/mosdepth_d4
- fraguracy: https://github.com/brentp/fraguracy/releases/download/v0.2.4/fraguracy

There is also a utility script `envs/download_containers_executables.sh` to download executables and [containers](#containers).

Add absolute paths to `snakemake_config.yaml`

### Containers 

Here is a list of containers required to run the Snakemake workflows

- bcftools: https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0
- bedtools_2.31.1: https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2
- bedtools: https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0
- d4tools: docker://clinicalgenomics/d4tools:2.0
- deepvariant: docker://quay.io/nf-core/deepvariant:1.5.0
- gatk: https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0
- happy: docker://community.wave.seqera.io/library/hap.py_rtg-tools:2ebb433f3ce976d3
- mosdepth: https://depot.galaxyproject.org/singularity/mosdepth:0.3.8--hd299d5a_0
- multiqc: https://depot.galaxyproject.org/singularity/multiqc:1.21--pyhdfd78af_0
- pandas: https://depot.galaxyproject.org/singularity/pandas:1.5.2
- picard: https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1
- rtgtools: docker://realtimegenomics/rtg-tools:3.12.1
- samtools: https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0

Pull images from the URLs or tag using singularity

```
singularity pull <image.img> <url/tag>
```

or Apptainer

```
apptainer pull <image.img> <url/tag>
```

There is also a utility script `envs/download_containers_executables.sh` to download executables and [containers](#containers), apptainer or singularity is specified using the `-t` flag, e.g. `bash download_containers_executables.sh -t apptainer`. The script downloads to the current directory. 

Finally add the absolute paths to the config `snakemake_config.yaml`

### Avidity Manuscript environment container

In the `envs/avidity` folder is a conda environment YAML and singularity definition file to generate a container required for running the following snakemake analysis workflows.

- `analysis/chr20_duplicates/Snakefile`
- `analysis/stack_reads/Snakefile`

To generate the required container using Singularity

```
cd envs/avidity
singularity build avidity.sif Singularity.def
```

or Apptainer

```
cd envs/avidity
apptainer build avidity.sif Singularity.def
```

Add the container image path to the `snakemake_config.yaml` config.


### Jupyther notebooks execution

Jupyter notebooks in the folder `notebooks` were executed using a conda as defined in `envs/environment.yml`. Create the conda environment using: 

```
conda env create -f envs/environment.yml
```

## Resources

### Genome download and setup

The reference genome can be downloaded and generate required files using a utility Snakemake workflow. `snakemake` v8.20.1 or later is required to run and the `samtools` container (see [Containers](#containers)) is needed for indexing. 

```
cd resources/GRCh38_GIABv3/
snakemake -j 4 -k -p --use-singularity
```

The FASTA must be un-bgzipped due to a sarek issue (see https://github.com/nf-core/sarek/issues/1741).


### GIAB stratifications

Genome in a bottle (GIAB) stratifications (v3.5) can be downloaded using the link below

https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/genome-stratifications-GRCh38@all.tar.gz

To download and extract the data in to the resource folder run:

```
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/genome-stratifications-GRCh38@all.tar.gz
tar -xzf genome-stratifications-GRCh38@all.tar.gz -C resources
```

### Avidity codebase

Some of the analysis relies on scripts developed for the paper [Arslan et al. Sequencing by avidity enables high accuracy with low reagent consumption. Nat Biotechnol 42, 132–138 (2024)](https://doi.org/10.1038/s41587-023-01750-7).

These are available at this repository: https://github.com/Elembio/AvidityManuscript2023 (specifically commit 701be395c892d00beca69693536ad600d209eec2).

Clone this into the `scripts` folder using the following command:

```
git clone --revision=701be395c892d00beca69693536ad600d209eec2 --depth=1 https://github.com/Elembio/AvidityManuscript2023.git scripts/AvidityManuscript2023
```

## Running workflows

### Running nf-core/sarek on Element/Illumina data

[Nextflow](http://nextflow.io) v24.10.4 is required for running sarek.

Download the [FASTQs](#illumina--element-fastq-download) and [genome reference](#genome-download).

The configs for running sarek are found in the `analysis/nfcore_sarek_rerun` folder with one subfolder for each data source. 

```
├── aviti_hq   # Element AVITI CB
├── aviti_ngi  # Element AVITI CB FS
└── xplus_sns  # Illumina NovaSeq XPLUS 
```

Within each is as `*params_relpath.yaml` defining the necessary parameters. Do not use the `*params.yaml` files as they are configured for a specific HPC resource. 

The `nextflow_no_tower.config` configures the FASTP step to cap reads to 150 bp but perform no other trimming. 

Modify the command below with your profile for choise (e.g. `singularity`) and run within the subfolder with the corresponding parameter YAML. E.g.

```
cd analysis/nfcore_sarek_rerun/xplus_sns
nextflow run nf-core/sarek -r 3.4.2 -profile <profile> -c ../nextflow_no_tower.config -params-file xplus_sns_params_relpath.yaml
```

### Running nf-core/sarek on public Element FreeStyle data

Besides the data generated to this study we also relied on public data, see e.g. `analysis/public_data`. One of the datasets require mapping, to download see [here](#public-element-freestyle-dataset). 

Mapping was performed using nf-core/sarek and parameters/configs/samplesheets are found in the `analysis/nfcore_sarek_public` folder. Follow instructions from [here](#running-nf-coresarek-on-elementillumina-data)

### Running secondary snakemake analysis

Snakefiles are found in the `analysis` folder and requires `snakemake` v8.20.1 or later to run. Most snakemake workflows requires [nf-core/sarek to be run first](#running-nf-coresarek-on-elementillumina-data). 

Configurations for resource files, executables and containers needed are specified in a YAML: `snakemake_config.yaml`. Make sure all paths are updated before you run.
