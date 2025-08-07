# Downloaded reference

https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/references/GRCh38/

# Build BWA index

sbatch -A ngi2016004 -J bwaindex -n 4 -p core -t 12:00:00 --wrap "apptainer exec /vulpes/ngi/containers/biocontainers/singularity-bwa-0.7.17--hed695b0_7.img bwa index genome.fa" --mail-user phojer@kth.se --mail-type ALL

# Unbgzipped index

Unzip

gunzip -c GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz > unbgzipped/genome.fa

Index FASTA

ml bioinfo-tools samtools
samtools faidx unbgzipped/genome.fa 

Build BWA index

cd unbgzipped
sbatch -A ngi2016004 -J bwaindex -n 4 -p core -t 12:00:00 --wrap "apptainer exec /vulpes/ngi/containers/biocontainers/singularity-bwa-0.7.17--hed695b0_7.img bwa index genome.fa" --mail-user phojer@kth.se --mail-type ALL


