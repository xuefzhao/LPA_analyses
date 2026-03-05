#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 input.tsv reference.fasta output.tsv"
    exit 1
fi

INPUT=$1
REF=$2
OUTPUT=$3

TMP_QUERY="tmp_query.fa"
TMP_BLAST="tmp_blast.tsv"

# -------------------------------------------------
# 1. Convert tab file (name \t sequence) to FASTA
# -------------------------------------------------
awk -F'\t' '{
    print ">"$1
    print $2
}' "$INPUT" > "$TMP_QUERY"

# -------------------------------------------------
# 2. Make BLAST database (if not already present)
# -------------------------------------------------
if [[ ! -f "${REF}.nin" && ! -f "${REF}.nhr" ]]; then
    ./ncbi-blast-2.17.0+/bin/makeblastdb -in "$REF" -dbtype nucl
fi

# -------------------------------------------------
# 3. Run BLAST
# -------------------------------------------------
./ncbi-blast-2.17.0+/bin/blastn \
    -task blastn-short \
    -query "$TMP_QUERY" \
    -db "$REF" \
    -outfmt "6 qseqid sseqid sstart send bitscore sstrand sseq" \
    -max_target_seqs 1000000 \
    -max_hsps 4 \
    > "$TMP_BLAST"

# -------------------------------------------------
# 4. Format final output
# -------------------------------------------------
echo -e "name\treference\tstart\tend\tbitscore\tsstrand\tsseq" > "$OUTPUT"

awk '$5 > 100 {
    start=$3
    end=$4
    if (start > end) {
        tmp=start
        start=end
        end=tmp
    }
    print $1"\t"$2"\t"start"\t"end"\t"$5"\t"$6"\t"$7
}' "$TMP_BLAST" >> "$OUTPUT"

# Cleanup
rm -f "$TMP_QUERY" "$TMP_BLAST"

echo "Done. Output written to $OUTPUT"





