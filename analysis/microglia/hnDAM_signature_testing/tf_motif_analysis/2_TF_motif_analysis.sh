# run in interactive session

cd /data/ADRD/glia_across_NDDs/analysis/microglia/DAM_signature_discovery/tf_motif_analysis

module load meme
module load samtools
module load bedtools

REF_DIR="refdata-gex-GRCh38-2024-A"
MOTIF_DIR="motifs"
OUT_DIR="outs"
PROMOTER_RANGE_UP=1000
PROMOTER_RANGE_DOWN=100



# === STEP 1: Extract TSS from compressed GTF ===
zcat "$REF_DIR/genes/genes.gtf.gz" | awk '$3 == "transcript"' | \
awk 'BEGIN{OFS="\t"} {
    split($10, gene_id, "\"");
    if ($7 == "+") {
        print $1, $4-1, $4, gene_id[2], ".", $7
    } else {
        print $1, $5-1, $5, gene_id[2], ".", $7
    }
}' > "$OUT_DIR/tss.bed"



# === STEP 2: Generate promoter regions ===
samtools faidx "$REF_DIR/fasta/genome.fa"
cut -f1,2 "$REF_DIR/fasta/genome.fa.fai" > "$OUT_DIR/genome.chrom.sizes"

bedtools flank -i "$OUT_DIR/tss.bed" \
               -g "$OUT_DIR/genome.chrom.sizes" \
               -l $PROMOTER_RANGE_UP -r $PROMOTER_RANGE_DOWN -s > "$OUT_DIR/promoters.bed"



# === STEP 3: Get FASTA sequences of promoter regions ===
bedtools getfasta -fi "$REF_DIR/fasta/genome.fa" \
                  -bed "$OUT_DIR/promoters.bed" -s -name > "$OUT_DIR/promoters.fa"




# === STEP 4: Run FIMO on each motif ===
for motif_file in "$MOTIF_DIR"/*.meme; do
    motif_name=$(basename "$motif_file" .meme)
    fimo_out="$OUT_DIR/fimo_${motif_name}"
    mkdir -p "$fimo_out"

    fimo --oc "$fimo_out" --thresh 1e-4 "$motif_file" "$OUT_DIR/promoters.fa"

    if [ -f "$fimo_out/fimo.tsv" ]; then
        echo "[*] Mapping motif hits back to genes for $motif_name..."
        awk 'BEGIN {OFS="\t"} NR > 1 && $3 ~ /^[0-9]+$/ && $4 ~ /^[0-9]+$/ { print $2, $3, $4, $1, $7, $5 }' "$fimo_out/fimo.tsv" > "$fimo_out/fimo.bed"
        bedtools intersect -a "$fimo_out/fimo.bed" -b "$OUT_DIR/promoters.bed" -wa -wb > "$fimo_out/${motif_name}_motif_hits.tsv"
    else
        echo "[!] Warning: No hits found for $motif_name"
    fi
done
