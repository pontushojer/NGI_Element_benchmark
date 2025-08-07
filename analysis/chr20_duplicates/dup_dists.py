import argparse
from collections import defaultdict

# https://gatk.broadinstitute.org/hc/en-us/articles/360036862931-MarkDuplicates-Picard
OPTICAL_DISTANCE = 2500 # recommended for patterned FCs

def main(args):
    dups = defaultdict(lambda: defaultdict(set))
    dups2 = defaultdict(set)

    dups_sizes = defaultdict(set)

    with open(args.sam, "r") as file:
        n = 0
        for line in file:
            if line.startswith("@"):
                continue
            
            fields = line.strip().split("\t")
            readname = fields[0]
            rg = fields[-3]
            di = fields[-2]
            ds = fields[-1]
            if n == 0:
                assert rg.startswith("RG:Z:")
                assert di.startswith("DI:i:")
                assert ds.startswith("DS:i:")
            
            rg = rg.split(":")[-1]
            di = int(di.split(":")[-1])
            ds = int(ds.split(":")[-1])

            *_, tile, x, y = readname.split(":")
            x = int(x)
            y = int(y)

            # Unique group
            group = (tile)

            if (x, y) not in dups[group][di]:
                n += 1
            dups[group][di].add((x, y))
            dups2[di].add((tile, x, y))
            dups_sizes[ds].add(di)    

    reldist = defaultdict(int)
    n_within_tile = 0
    n_optical = 0 
    for group, di in dups.items():
        for di, coords in di.items():
            coords = list(coords)

            # Skip if only one coordinate in the group meaning that they 
            # are not on the same tile
            if len(coords) == 1:
                continue

            n_within_tile += len(coords)
            #assert len(coords) > 1, f"Only one coordinate for group {group} and DI {di}: {coords}"
            while len(coords) > 1:
                dists = []
                x1, y1 = coords.pop()
                for x2, y2 in coords:
                    # Calculate euclidean distance
                    # This is not the same condition used in MarkDuplicates which simply checks if 
                    # the distance in x and y is less than the threshold (default 100)
                    # see: https://github.com/broadinstitute/picard/blob/c8523a0dfbe097a91ed8ffb0b0972cb7f9395595/src/main/java/picard/sam/markduplicates/util/OpticalDuplicateFinder.java#L337
                    dist = ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5
                    dists.append(dist)

                min_dist = min(dists)
                if min_dist < OPTICAL_DISTANCE:
                    n_optical += 1
                reldist[100 * (int(min_dist)//100)] += 1
    
    for di, reads in dups2.items():
        ds = len(reads)
        if di not in dups_sizes[ds]:
            for size in dups_sizes:
                if di in dups_sizes[size]:
                    print(di, size, ds)

                


    print(f"Number of reads: {n}")
    print(f"Number of reads within tile: {n_within_tile} ({n_within_tile / n:.2%})")
    print(f"Number of optical (dist<{OPTICAL_DISTANCE}): {n_optical} ({n_optical / n:.2%})")
    
    print(sum(reldist.values()))
    print("\nDistance\tCount\tCumulativeFraction")
    total = 0
    for dist, count in sorted(reldist.items()):
        total += count
        frac_within_distance = total / n_within_tile
        print(f"{dist}\t{count}\t{frac_within_distance:.1%}")


    total = sum(map(len, dups_sizes.values()))
    print(total)
    print("\nSize\tCount\tPercentGroups\tPercentReads")
    for size, dis in sorted(dups_sizes.items()):
        count = len(dis)
        fraction = count / total
        frac_reads = count * size / n
        print(f"{size}\t{count}\t{fraction:.1%}\t{frac_reads:.1%}")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze duplicate distances.")
    parser.add_argument("sam", type=str, help="Input SAM path with only RG:Z, DI:i and DS:i tags")
    args = parser.parse_args()
    
    main(args)
