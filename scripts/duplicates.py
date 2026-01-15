import argparse
from collections import defaultdict, Counter
from itertools import combinations
from dataclasses import dataclass
from math import hypot

import pysam 

@dataclass 
class Position:
    tile: str
    x: int
    y: int

    def distance(self, other):
        if self.tile != other.tile:
            return float("inf")
        dx = abs(self.x - other.x)
        dy = abs(self.y - other.y)
        return hypot(dx, dy)

@dataclass(eq=False)
class ReadPair:
    readf: pysam.AlignedSegment
    readr: pysam.AlignedSegment

    @classmethod
    def from_reads(cls, reads):
        assert len(reads) == 2
        
        # Assign by relative order
        if not reads[0].is_reverse and reads[1].is_reverse:
            return cls(reads[0], reads[1])
        elif reads[0].is_reverse and not reads[1].is_reverse: 
            return cls(reads[1], reads[0])        
        else:
            raise ValueError("Not inward facing")

    def get_position(self):
        name = self.readf.query_name
        *_, tile, x, y = name.split(":")
        x = int(x)
        y = int(y)        
        return Position(tile, x, y)

    def is_fr(self):
        return self.readf.is_read1

    @property
    def name(self):
        return self.readf.query_name
    
    def __eq__(self, other):
        return self.name == other.name
    
    def __hash__(self):
        return hash(self.name)

    def distance_to(self, other):
        pos1 = self.get_position()
        pos2 = other.get_position()
        return pos1.distance(pos2)

def main(args):
    stats = Counter()

    dups = defaultdict(list)
    with pysam.AlignmentFile(args.bam, 'rb') as file:
        for read in file.fetch():
            # Only duplicates with DI tag
            if not read.has_tag("DI"):
                continue

            if read.is_supplementary:
                continue

            # Skip non paired reads
            if not read.is_proper_pair:
                continue
            
            stats["Reads"] += 1

            di = read.get_tag("DI")
            
            dups[di].append(read)
    
    dups_by_dist = defaultdict(Counter)
    for di, reads in dups.items():
        # Group by read name to get pairs
        reads_by_name = defaultdict(list)
        for read in reads:
            reads_by_name[read.query_name].append(read)
        
        # Create list of ReadPair objects
        pairs = []
        for reads in reads_by_name.values():
            pair = ReadPair.from_reads(reads)
            pairs.append(pair)
            stats["Read pairs"] += 1

        if len(pairs) == 2:
            stats["Read pairs ds=2"] += 2

        # Find pairs that are optical duplicates
        optical = set()
        for p1, p2 in combinations(pairs, 2):
            stats["Pair combos"] += 1

            distance = p1.distance_to(p2)
            proximal = distance < args.dist_optical
            complementary = p1.is_fr() != p2.is_fr()

            if proximal and not complementary:
                optical |= {p1, p2}
            
            ptype = "Prox" if proximal else "noProx"
            ptype += "+Compl" if complementary else "+noCompl"

            stats[f"Pair combos {ptype}"] += 1

            if distance != float("inf"):
                distance_by_10 = 10 * (distance // 10)
                dups_by_dist[distance_by_10][complementary] += 1 

    for k in list(stats):
        if k.startswith("Pair combos "):
            stats["Pct" + k] = round(100 * stats[k] / stats["Pair combos"], 1)
    
    with open(args.prefix + ".stats.txt", "w") as f:
        for k,v in sorted(stats.items()):
            print(f"{k:<30}: {v:>10}", file=f)


    with open(args.prefix + ".distances.csv", "w") as f:
        print("Distance,Dups,ComplementaryDups,NoncomplementaryDups", file=f)
        for dist, data in sorted(dups_by_dist.items()):
            comp_dups = data[True]
            noncomp_dups = data[False]
            total_dups = comp_dups + noncomp_dups

            print(dist,total_dups,comp_dups,noncomp_dups, sep=",", file=f)
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze duplicate types.")
    parser.add_argument("bam", type=str, help="Input BAM path with DI:i tags from MarDuplicates.")
    parser.add_argument("-d", "--dist-optical", default=2500,
                        help="Physical distance above which reads are consider optical duplicates. "
                             "See https://gatk.broadinstitute.org/hc/en-us/articles/360036862931-MarkDuplicates-Picard. "
                              "Default: %(default)s (recommended for patterned FCs)")
    parser.add_argument("-p", "--prefix", 
                        help="Prefix for output files ´<prefix>.stats.txt´ and ´<prefix>.distances.csv´")
    args = parser.parse_args()
    
    main(args)
