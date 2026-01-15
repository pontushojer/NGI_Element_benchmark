# G4 data

## pqsfinder data

https://pqsfinder.fi.muni.cz/genomes

All tracks were computed using pqsfinder v2.0.1 with default algorithm options.

- https://pqsfinder.fi.muni.cz/hub/hg38/pqsfinder_hg38.bb
- https://pqsfinder.fi.muni.cz/hub/hg38/pqsfinder_hg38.gff
- https://pqsfinder.fi.muni.cz/hub/hg38/pqsfinder_hg38.bed

Sorted with `bedtools sort`

Converted to 3 field BED and bgzipped using:

cut -f -3 pqsfinder_hg38.sort.bed | bgzip -c > pqsfinder_hg38.sort.3tab.bed.gz