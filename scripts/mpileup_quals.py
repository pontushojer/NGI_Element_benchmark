"""
Script to extract quality scores from a mpileup file.

For each continous region the quality scores histogram for the relative position in the region.
"""

import argparse
from collections import Counter, defaultdict


def ascii_to_phred(ascii):
    return ord(ascii) - 33


def parse_bed(bed_file):
    regions = []
    with open(bed_file) as f:
        for line in f:
            els = line.strip().split('\t')
            if len(els) == 3:
                chrom, start, end = els
                strand = "+"
            else:
                chrom, start, end, _, _, strand = els

            regions.append((chrom, int(start), int(end), strand))
    return regions

import re


def extract_bases(pileup_str):
    """
    Extract the actual bases covering a pileup position from a pileup string.
    
    Args:
        pileup_str (str): The pileup base string from tools like samtools mpileup.
        
    Returns:
        list: A list of base characters (A, C, G, T, ., ,) observed at the position.
              '.' and ',' represent matches to the reference on forward and reverse strands.
    """
    bases = []
    i = 0
    length = len(pileup_str)
    
    while i < length:
        char = pileup_str[i]

        if char == '^':
            # Start of a read segment, skip the next character (mapping quality)
            i += 2
        elif char == '$':
            # End of read segment
            i += 1
        elif char in '+-':
            # Indel detected; skip over the length of the indel and the bases
            bases.append(char)
            i += 1
            indel_len_str = ''
            while i < length and pileup_str[i].isdigit():
                indel_len_str += pileup_str[i]
                i += 1
            indel_len = int(indel_len_str)
            i += indel_len  # Skip over the indel sequence
        else:
            # Actual base
            bases.append(char.upper())
            i += 1

    return bases

assert extract_bases("TTT.TT..T..").count("T") == 6, extract_bases("TTT.TT..T..")
assert len(extract_bases("TTT.TT..T..")) == 11, extract_bases("TTT.TT..T..")
assert extract_bases(".+5CCATG.+5CCATG.+5CCATG..+5CCATG.+5CCATG..G..").count("+") == 5, extract_bases(".+5CCATG.+5CCATG.+5CCATG..+5CCATG.+5CCATG..G..")
assert extract_bases(".+5CCATG.+5CCATG.+5CCATG..+5CCATG.+5CCATG..G..").count("G") == 1, extract_bases(".+5CCATG.+5CCATG.+5CCATG..+5CCATG.+5CCATG..G..")
assert len(extract_bases(".+5CCATG.+5CCATG.+5CCATG..+5CCATG.+5CCATG..G..")) == 16, extract_bases(".+5CCATG.+5CCATG.+5CCATG..+5CCATG.+5CCATG..G..")

def parse_pileup_region(pileup_file):
    with open(pileup_file) as f:
        region = []
        
        last_postion = None
        for line in f:
            chrom, position, refbase, n_reads, bases, quals = line.strip().split('\t')
            n_reads = int(n_reads)
            position = int(position)

            if last_postion is not None and position != (last_postion + 1):
                yield region
                region = []

            # No reads
            if n_reads == 0:
                quals = []
                bases = []
            else:
                quals = list(quals)
                bases = extract_bases(bases)
            
            assert len(quals) == n_reads <= len(bases), (len(quals), n_reads, len(bases), line)

            region.append((chrom, position, quals, bases, int(n_reads)))
            
            last_postion = position

        yield region    

def compute_phread_histogram(hist):
    hist_phred = defaultdict(Counter)
    
    max_phred = 0
    for pos, counter in hist.items():
        for qual, count in counter.items():
            phred = ascii_to_phred(qual)
            max_phred = max(max_phred, phred)
            hist_phred[pos][phred] = count

    return hist_phred, max_phred


def print_histogram(hist, is_same_strand, phread_span, nreads_pos):
    for pos, phread_hist in sorted(hist.items()):
        n_reads = nreads_pos[pos]
        pos_quals = [phread_hist[phred] for phred in phread_span]
        print(pos, is_same_strand, n_reads, *pos_quals, sep='\t')


def print_bases_histogram(hist, is_same_strand, bases, nreads_pos):
    for pos, base_hist in sorted(hist.items()):
        n_reads = nreads_pos[pos]
        #assert n_reads <= sum(base_hist.values()), (n_reads, sum(base_hist.values()))
        pos_bases = [base_hist[base] for base in bases]
        print(pos, is_same_strand, n_reads, *pos_bases, sep='\t')


def relative_position_middle(position, start, end, reverse):
    """Center on the feature and report relative offset"""
    # Inside feature
    middle = (start + end) // 2

    # Calculate the distance from the middle of the feature
    if reverse:
        return middle - position
    return position - middle


def relative_position_squach(position, start, end, reverse):
    """Squach all positions inside feature to 0 and report rest relative to that"""
    if position < start:
        dist = abs(start - position)
        if reverse:
            return dist
        return -dist
    
    elif position > end:
        dist = abs(end - position)
        if reverse:
            return -dist
        return dist

    else: # Inside feature
        return 0

def relative_position_restrict(position, start, end, reverse, linside=10):
    """Squach all positions inside feature to a `linside` window report rest relative to that"""
    assert linside % 2 == 0

    # Used of offset positions
    half = linside // 2

    if position < start:
        dist = abs(start - position)
        if reverse:
            return dist + half
        return -dist - half
    
    elif position > end:
        dist = abs(end - position)
        if reverse:
            return -dist - half
        return dist + half

    else: # Inside feature
        dist = end - position if reverse else position - start
        l = end-start
        return int(linside * dist / l - half)


def get_relative_position(position, start, end, reverse, method=None):
    match method:
        case "middle" | None: # Default
            return relative_position_middle(position, start, end, reverse)
        case "squach":
            return relative_position_squach(position, start, end, reverse)
        case "restrict":
            return relative_position_restrict(position, start, end, reverse)
        case _:
            raise ValueError(f"Unknown method {method}")


def main(args):
    hist_ss = defaultdict(Counter)
    hist_ds = defaultdict(Counter)
    nreads_ds = Counter()
    nreads_ss = Counter()
    for pileup_region, feature in zip(parse_pileup_region(args.pileup), parse_bed(args.bed)):
        chrom, start, end, strand = feature

        if not pileup_region:
            continue
            
        try:
            assert chrom == pileup_region[0][0]
            assert pileup_region[0][1] < start
            assert end < pileup_region[-1][1]
        except AssertionError as e:
            print(e)
            print(feature)
            print(pileup_region[0][0], pileup_region[0][1], pileup_region[-1][1])
            exit(1)


        feature_reverse = strand == '-'
        same_strand = feature_reverse == args.reverse

        for chrom, position, quals, bases, nreads in pileup_region:
            relative_pos = get_relative_position(position, start, end, args.reverse, method=args.distance_method)
            
            if same_strand:
                nreads_ss[relative_pos] += nreads
            else:
                nreads_ds[relative_pos] += nreads

            if not args.mismatch_bases:
                if same_strand:
                    hist_ss[relative_pos].update(quals)
                else:
                    hist_ds[relative_pos].update(quals)
            else:
                if same_strand:
                    hist_ss[relative_pos].update(bases)
                else:
                    hist_ds[relative_pos].update(bases)
    
        
    if not args.mismatch_bases:
        hist_ds, max_phred_ds = compute_phread_histogram(hist_ds)
        hist_ss, max_phred_ss = compute_phread_histogram(hist_ss)

        max_phred = max(max_phred_ds, max_phred_ss)
        phread_span = list(range(max_phred + 1))
        # Print header with
        # - relative position =0 is inside the feature, <0 is upstream, >0 is downstream
        # - ss - same strand as the feature
        # - number of bases
        # - quality scores from phred 0 to max_phred

        print('pos', 'ss', 'bases', *phread_span, sep='\t')
        print_histogram(hist_ds, False, phread_span, nreads_ds)
        print_histogram(hist_ss, True, phread_span, nreads_ss)
        

    else:
        bases = list('ACGT+-')
        print('pos', 'ss', 'bases', *bases, sep='\t')
        print_bases_histogram(hist_ds, False, bases, nreads_ds)
        print_bases_histogram(hist_ss, True, bases, nreads_ds)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pileup', help='BED file with pileups')
    parser.add_argument('bed', help='BED file with regions')
    parser.add_argument('-r', '--reverse', action='store_true', help='Reads are on the reverse strand')
    parser.add_argument("-f", "--flank", type=int, default=150, help="Flank size")
    parser.add_argument("-m", "--mismatch-bases", action='store_true', help="Output mismatch bases instead of quality scores")
    parser.add_argument("-d", "--distance-method", choices=["middle", "squach", "restrict"], default="middle",
                        help="Method for reporting distances to feature. " \
                            "'middle' = center on feature (default), " \
                            "'squach' = positions inside feature reported as 0, " \
                            "'restrict' = positions inside feature limited to 10 base from which the outside positions others are offset.")
    args = parser.parse_args()
    main(args)