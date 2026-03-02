#!/usr/bin/env python3

import sys
from collections import defaultdict

if len(sys.argv) != 5:
    print("Usage: python extract_spanned_regions.py reference.fa blast_output.tsv flank_size output.fa")
    sys.exit(1)

fasta_file = sys.argv[1]
blast_file = sys.argv[2]
flank_size = int(sys.argv[3])
output_file = sys.argv[4]

# ----------------------------------------
# 1. Parse BLAST output → get min/max per reference
# ----------------------------------------
regions = defaultdict(lambda: [None, None])

with open(blast_file) as f:
    header = f.readline()  # skip header
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) < 4:
            continue

        ref = fields[1]
        start = int(fields[2])
        end = int(fields[3])

        if regions[ref][0] is None or start < regions[ref][0]:
            regions[ref][0] = start
        if regions[ref][1] is None or end > regions[ref][1]:
            regions[ref][1] = end

# ----------------------------------------
# 2. FASTA reader
# ----------------------------------------
def fasta_reader(fp):
    header = None
    seq = []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if header:
                yield header, "".join(seq)
            header = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)
    if header:
        yield header, "".join(seq)

# ----------------------------------------
# 3. Extract sequences with flanking
# ----------------------------------------
with open(output_file, "w") as out:
    with open(fasta_file) as f:
        for header, sequence in fasta_reader(f):
            seq_len = len(sequence)

            if header in regions and regions[header][0] is not None:
                start, end = regions[header]

                # Add flanking
                new_start = max(1, start - flank_size)
                new_end = min(seq_len, end + flank_size)

                # Convert to 0-based indexing
                subseq = sequence[new_start - 1:new_end]

                out.write(f">{header}:{new_start}-{new_end}\n")
                out.write(subseq + "\n")
            else:
                # No hits → keep entire sequence
                out.write(f">{header}\n")
                out.write(sequence + "\n")

print("Done.")