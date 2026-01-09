[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17302040.svg)](https://doi.org/10.5281/zenodo.17302040)
# Element vs Illumina benchmark

Project: SR_23_02

## Data organisation

- data/wgs: Raw FASTQs. To download data, see [Data download](#data-download) section.
- analysis: Analysis runs on data, either using Snakemake script or existing Nextflow/nf-core workflow
- resources: Genome annotations and reference
- notebooks: Jupyther notebooks for data processing and visualization
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
