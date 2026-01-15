#!/bin/bash
# Download containers and executables into the current directory

# Default value
tool="singularity"

while getopts "t:" opt; do
  case $opt in
    t)
      # Check if the argument is valid
      if [[ "$OPTARG" == "singularity" || "$OPTARG" == "apptainer" ]]; then
        tool="$OPTARG"
      else
        echo "Error: Invalid type '$OPTARG'. Must be 'singularity' or 'apptainer'." >&2
        exit 1
      fi
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

echo "Configured to use: $tool"

set -euox pipefail

wget -O mosdepth_d4_0.3.10 https://github.com/brentp/mosdepth/releases/download/v0.3.10/mosdepth_d4
chmod +x mosdepth_d4_0.3.10

wget -O fraguracy_0.2.4 https://github.com/brentp/fraguracy/releases/download/v0.2.4/fraguracy
chmod +x fraguracy_0.2.4

$tool pull singularity-bcftools-1.18--h8b25389_0.img https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0
$tool pull singularity-bedtools-2.31.1--hf5e1c6e_2.img https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2
$tool pull singularity-bedtools-2.30.0--hc088bd4_0.img https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0
$tool pull clinicalgenomics-d4tools.sif docker://clinicalgenomics/d4tools:2.0
$tool pull nf-core-deepvariant-1.5.0.img docker://quay.io/nf-core/deepvariant:1.5.0
$tool pull singularity-gatk4-4.5.0.0--py36hdfd78af_0.img https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0
$tool pull hap.py_rtg-tools_2ebb433f3ce976d3.img docker://community.wave.seqera.io/library/hap.py_rtg-tools:2ebb433f3ce976d3
$tool pull singularity-mosdepth-0.3.8--hd299d5a_0.img https://depot.galaxyproject.org/singularity/mosdepth:0.3.8--hd299d5a_0
$tool pull singularity-multiqc-1.21--pyhdfd78af_0.img https://depot.galaxyproject.org/singularity/multiqc:1.21--pyhdfd78af_0
$tool pull depot.galaxyproject.org-singularity-pandas-1.5.2.img https://depot.galaxyproject.org/singularity/pandas:1.5.2
$tool pull singularity-picard-3.0.0--hdfd78af_1.img https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1
$tool pull rtg-tools_realtimegenomics_3.12.1.sif docker://realtimegenomics/rtg-tools:3.12.1
$tool pull singularity-samtools-1.19.2--h50ea8bc_0.img https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0
