[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17302040.svg)](https://doi.org/10.5281/zenodo.17302040)
# Element vs Illumina benchmark

Project: SR_23_02

## Data organisation

- data/wgs: Raw FASTQs organized under here. To download data, see [Data download](#data-download) section.
- analysis: Analysis runs on data, either using Snakemake script or existing Nextflow/nf-core workflow
- resources: Genome annotations and reference
- notebooks: Jupyter notebooks for data processing and visualization
- scripts: Custom Python scripts
- figures/svg: SVG figures

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

## Genome download 

The reference genome can be downloaded using a utility script

```
cd resources/GRCh38_GIABv3/
bash download_and_unzip.sh
```

The FASTA is un-bgzipped due to a sarek issue (see https://github.com/nf-core/sarek/issues/1741).


## Running nf-core/sarek on Element/Illumina data

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


## Jupyther notebooks execution

Jupyter notebooks in the folder `notebooks` were executed using a conda as defined in `environment.yml`. Create the conda environment using: 

```
conda env create -f environment.yml
```


