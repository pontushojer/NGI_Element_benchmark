"""
Script to extract quality scores from a mpileup file.

For each continous region the quality scores histogram for the relative position in the region.
"""
import argparse
from collections import defaultdict

from mpileup_quals import parse_pileup_region, parse_bed, get_relative_position, ascii_to_phred


def main(args):
    stats_per_region = {
        True: defaultdict(lambda: defaultdict(dict)), # same strand
        False: defaultdict(lambda: defaultdict(dict)) # different strand
    }
    feature_coverage = defaultdict(lambda: defaultdict(dict))

    features = list(parse_bed(args.bed))
    pileup_names = [pileup_file.split('/')[-1] for pileup_file in args.pileups]
    assert len(set(pileup_names)) == len(pileup_names)

    for pileup_file, pileup_name in zip(args.pileups, pileup_names):
        print(f'Processing {pileup_file}')
        n = 0
        for pileup_region, feature in zip(parse_pileup_region(pileup_file), features):
            chrom, start, end, strand = feature
            feature_name = f"{chrom}:{start}-{end}({strand})"
            if not pileup_region:
                continue
        
            assert chrom == pileup_region[0][0]
            
            feature_reverse = strand == '-'
            same_strand = feature_reverse == args.reverse
            stats = {}
            max_coverage = 0
            for chrom, position, quals in pileup_region:
                max_coverage = max(max_coverage, len(quals))
                if not quals:
                    continue

                relative_pos = get_relative_position(position, start, end, args.reverse)
                
                qs = [ascii_to_phred(q) for q in quals]
                
                # Calculate the percentage of quality scores greater than 30
                gt_q30 = sum(1 for q in qs if q >= 30)
                pct_gt_q30 = 100 * gt_q30 / len(qs)
                stats[relative_pos] = round(pct_gt_q30, 1)
                
                # Median quality score - not good for Illumina Q-binned data
                # median_q = sorted(qs)[len(qs) // 2]
                # stats[relative_pos] = median_q   

            feature_coverage[feature_name][pileup_name] = max_coverage
            stats_per_region[same_strand][feature_name][pileup_name] = stats
            n += 1
        print(f'Processed {n} regions')
    
    print('Finished processing pileup files')

    with open(args.output, "w") as out:
        print('region', 'same_strand', 'relative_position', *pileup_names, sep='\t', file=out)
        for same_strand, regions in stats_per_region.items():
            for region_name, pileup_stats in regions.items():
                coverage = feature_coverage[region_name]
                if min(coverage.values()) <= args.min_coverage:
                    print(f"Skipping {region_name} due to low coverage")
                    continue

                # Reorganize stats to list the pileup with stats for each base
                stats_per_base = defaultdict(dict)
                for pileup_file, stats in pileup_stats.items():
                    for i in range(-args.flank, args.flank + 1):
                        if i not in stats:
                            continue
        
                        stats_per_base[i][pileup_file] = stats[i]

                min_pos = min(stats_per_base.keys())
                max_pos = max(stats_per_base.keys())
                for i in range(min_pos, max_pos + 1):
                    file_stats = stats_per_base.get(i, {})
                    pileup_stats = [file_stats.get(pileup_file, float('nan')) for pileup_file in pileup_names]
                    
                    print(region_name, same_strand, i, *pileup_stats, sep='\t', file=out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('bed', help='BED file with regions')
    parser.add_argument('pileups', help='BED files with pileups', nargs='+')
    parser.add_argument('-r', '--reverse', action='store_true', help='Reads are on the reverse strand')
    parser.add_argument("-f", "--flank", type=int, default=150, help="Flank size")
    parser.add_argument('-o', '--output', default='mpileup_quals_region.tsv', help='Output file')
    parser.add_argument('-m', '--min-coverage', type=int, default=10, help='Minimum coverage to include region')
    args = parser.parse_args()
    main(args)