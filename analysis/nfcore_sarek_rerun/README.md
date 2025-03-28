
# nf-core/sarek pipeline 2nd run

## Changes

Run CNV calling

```
-tools: deepvariant
+tools: deepvariant,cnvkit
```

Use different GRCh38 reference genome without ALTs and that masks some false duplications, source https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/
See also:
- https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02863-7 
- https://x.com/brent_p/status/1849503122114347051
- https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use


```
-genome: GATK.GRCh38
-igenomes_ignore: false
+fasta: /proj/ngi2016004/private/strategic_proj/SR_23_02_Element_vs_Illumina/resources/GRCh38_GIABv3/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
+fasta_fai: /proj/ngi2016004/private/strategic_proj/SR_23_02_Element_vs_Illumina/resources/GRCh38_GIABv3/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz.fai
+igenomes_ignore: true
+bwa: "/proj/ngi2016004/private/strategic_proj/SR_23_02_Element_vs_Illumina/resources/GRCh38_GIABv3/BWAIndex"
```

Skip base recalibration as its is not recommended for DeepVariant

```
+skip_tools: baserecalibrator,fastqc
```

## Run commands
Run in screen

Setup env
```
source /vulpes/ngi/production/latest/conf/sourceme_sthlm.sh
source activate NGI
```

Run nextflow

```
nextflow run /vulpes/ngi/production/v24.07/sw/sarek/3_4_2/ -profile uppmax --project ngi2016004 -c /vulpes/ngi/production/v24.07/conf/sarek_sthlm.config -c ../nextflow.config -params-file <params.file>
```

Use params.file from below

aviti_hq_params.yaml
aviti_ngi_params.yaml
xplus_sns_params.yaml



