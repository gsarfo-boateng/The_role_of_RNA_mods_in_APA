#!/bin/bash
set -euo pipefail

# Reference GTF and tools
REF_GTF="/mnt/george_drive/george/m6a/dm6_new.gtf"
ST_BIN="/mnt/george_drive/george/direct_RNAseq/stringtie-3.0.1.Linux_x86_64/stringtie"
SUPPA="/mnt/george_drive/george/direct_RNAseq/SUPPA-2.4/suppa.py"
THREADS=32

# Input BAMs
# NOTE: sample IDs are derived with: basename "$BAM" .minimap2.sorted.bam
# bcos the files end with ".sorted.bam", the SAMPLE IDs will include ".sorted.bam".
BAMS=(
"/mnt/george_drive/george/direct_RNAseq/W1118_induro_mRNA_2025_102/W1118_induro_mRNA_2025_102/20251002_1328_MN29498_FBE06035_a9b4fc34/pod5/W1118_induro_mRNA_2025_102.sorted.bam"
"/mnt/george_drive/george/direct_RNAseq/w1118_ind_rev_rep1_mRNA/w1118_ind_rev_rep1_mRNA/20250821_1257_MN29498_FBC55934_ad449135/pod5/w1118_ind_rev_rep1_mRNA_m6A_m5C_psi.sorted.bam"

"/mnt/george_drive/george/direct_RNAseq/FR113N_rep2/FR113N_rep2/FR113N_rep2/20250507_1202_MN29498_FBB33215_0fe8b5ba/pod5_skip/FR113N_05_11_2025_rep2_ref.sorted.bam"
"/mnt/george_drive/george/direct_RNAseq/FR113N_rep3/FR113N_rep3/20250508_1211_MN29498_FBB33215_84ccd776/pod5_skip/FR113N_05_11_2025_rep3_ref.sorted.bam"

#"/mnt/george_drive/george/direct_RNAseq/Mettl3_KO_mRNA_004_rep2/Mettl3_KO_mRNA_004_rep2/20250912_1021_MN29498_FBC38312_0de9666a/pod5/Mettl3_KO_mRNA_rep2.sorted.bam"
#"/mnt/george_drive/george/direct_RNAseq/Mettl3_KO_mRNA_RNA004/Mettl3_KO_mRNA_RNA004/20250910_1505_MN29498_FBC38311_72234379/pod5/Mettl3_KO_mRNA_RNA004.sorted.bam"
)

# =========================================
# CHOOSE ONE COMPARISON PER RUN
# (Run the script once per comparison)
# =========================================

# --- Run 1: KO vs W1118 ---
#COMP_NAME="KO_vs_W1118"
#FR_SAMPLES=("Mettl3_KO_mRNA_rep2.sorted.bam" "Mettl3_KO_mRNA_RNA004.sorted.bam")  # KO
#W_SAMPLES=("W1118_induro_mRNA_2025_102.sorted.bam" "w1118_ind_rev_rep1_mRNA_m6A_m5C_psi.sorted.bam")

# --- Run 2: KO vs FR113N (to use this, comment the Run 1 block above and uncomment below) ---
# COMP_NAME="KO_vs_FR113N"
# FR_SAMPLES=("Mettl3_KO_mRNA_rep2.sorted.bam" "Mettl3_KO_mRNA_RNA004.sorted.bam")  # KO
# W_SAMPLES=("FR113N_05_11_2025_rep2_ref.sorted.bam" "FR113N_05_11_2025_rep3_ref.sorted.bam")


# --- Run 3: FR113N vs w1118 (to use this, comment the Run 1 block above and uncomment below) ---
 COMP_NAME="FR113N_vs_w1118"
 FR_SAMPLES=("FR113N_05_11_2025_rep2_ref.sorted.bam" "FR113N_05_11_2025_rep3_ref.sorted.bam")  # KO
 W_SAMPLES=("W1118_induro_mRNA_2025_102.sorted.bam" "w1118_ind_rev_rep1_mRNA_m6A_m5C_psi.sorted.bam")


# Event types to analyze
EVENTS=("SE" "A3" "A5" "MX" "RI" "AF" "AL")

# =========================================
# PREP
# =========================================

mkdir -p as_work/{stringtie_asm,stringtie_quant,tpm,AS_event,psi,logs}
# diff dir now includes comparison name
mkdir -p "as_work/diff/${COMP_NAME}"

BASE_DIR="$(pwd -P)"

abspath () {
  if command -v realpath >/dev/null 2>&1; then
    realpath -m "$1"
  else
    python3 - <<'PY' "$1"
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
  fi
}

LOGDIR="$BASE_DIR/as_work/logs"
ASM_DIR="$BASE_DIR/as_work/stringtie_asm"
QNT_DIR="$BASE_DIR/as_work/stringtie_quant"
TPM_DIR="$BASE_DIR/as_work/tpm"
AS_EVT_DIR="$BASE_DIR/as_work/AS_event"
PSI_DIR="$BASE_DIR/as_work/psi"
DIFF_DIR="$BASE_DIR/as_work/diff/${COMP_NAME}"

command -v "$ST_BIN" >/dev/null 2>&1 || { echo "ERROR: stringtie not found at $ST_BIN"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found"; exit 1; }
command -v gawk  >/dev/null 2>&1 || { echo "ERROR: gawk required"; exit 1; }

# =========================================
# 1) StringTie per-sample assembly (long-read mode)
# =========================================
echo "[1/8] StringTie per-sample assembly..."
for BAM in "${BAMS[@]}"; do
  [ -f "$BAM" ] || { echo "Missing BAM: $BAM"; exit 1; }
  SAMPLE=$(basename "$BAM" .minimap2.sorted.bam)
  echo "  assembling: $SAMPLE"
  "$ST_BIN" -L -G "$REF_GTF" -p "$THREADS" \
    -o "$ASM_DIR/${SAMPLE}.gtf" "$BAM" \
    >"$LOGDIR/${SAMPLE}.stringtie_asm.log" 2>&1
done

# =========================================
# 2) StringTie merge
# =========================================
echo "[2/8] StringTie merge..."
ls "$ASM_DIR"/*.gtf > "$ASM_DIR/mergelist.txt"
"$ST_BIN" --merge -G "$REF_GTF" -p "$THREADS" \
  -o "$ASM_DIR/merged.gtf" "$ASM_DIR/mergelist.txt" \
  >"$LOGDIR/stringtie_merge.log" 2>&1
MERGED="$ASM_DIR/merged.gtf"

# =========================================
# 3) StringTie quant per sample vs merged.gtf (NO Ballgown)
# =========================================
echo "[3/8] StringTie quantification..."
for BAM in "${BAMS[@]}"; do
  SAMPLE=$(basename "$BAM" .minimap2.sorted.bam)
  echo "  quantifying: $SAMPLE"
  "$ST_BIN" -L -e -G "$MERGED" -p "$THREADS" \
    -A "$QNT_DIR/${SAMPLE}.gene_abund.tsv" \
    -o "$QNT_DIR/${SAMPLE}.merged.gtf" \
    "$BAM" >"$LOGDIR/${SAMPLE}.stringtie_quant.log" 2>&1
done

# =========================================
# 4) Build SUPPA-style TPM matrix from *.merged.gtf (Python)
# =========================================
echo "[4/8] Building SUPPA-style TPM matrix (Python)..."
TPM_OUT="$TPM_DIR/tpm_matrix.suppa.tsv"

python3 - "$QNT_DIR" "$TPM_OUT" <<'PY'
import sys, os, re, glob
qnt_dir = sys.argv[1]
out_path = sys.argv[2]
gtfs = sorted(glob.glob(os.path.join(qnt_dir, "*.merged.gtf")))
if not gtfs:
    sys.exit("No .merged.gtf files found in " + qnt_dir)

def parse_attrs(attr):
    def grab(k):
        m = re.search(rf'{k}\s+"([^"]+)"', attr)
        return m.group(1) if m else None
    return grab("transcript_id"), grab("TPM"), grab("FPKM")

sample_names = []
per_sample = []
all_tids = set()

for gtf in gtfs:
    sname = re.sub(r"\.merged\.gtf$", "", os.path.basename(gtf))
    sample_names.append(sname)
    tpm_by_tid = {}
    fpkm_by_tid = {}
    fpkm_sum = 0.0
    with open(gtf, "r", encoding="utf-8", errors="ignore") as fh:
        for ln in fh:
            if not ln or ln[0] == "#": continue
            cols = ln.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript": continue
            tid, tpm, fpkm = parse_attrs(cols[8])
            if not tid: continue
            if tpm is not None:
                try: tpm_by_tid[tid] = float(tpm)
                except: pass
            elif fpkm is not None:
                try:
                    v = float(fpkm)
                    fpkm_by_tid[tid] = v
                    fpkm_sum += v
                except: pass
    if not tpm_by_tid and fpkm_by_tid and fpkm_sum > 0:
        scale = 1_000_000.0 / fpkm_sum
        tpm_by_tid = {tid: v * scale for tid, v in fpkm_by_tid.items()}
    per_sample.append(tpm_by_tid)
    all_tids.update(tpm_by_tid.keys())

os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
with open(out_path, "w", encoding="utf-8") as out:
    out.write("\t".join(sample_names) + "\n")
    for tid in sorted(all_tids):
        row = [tid] + [str(d.get(tid, 0.0)) for d in per_sample]
        out.write("\t".join(row) + "\n")

print(f"Wrote {out_path} with {len(sample_names)} samples and {len(all_tids)} transcripts.")
PY

echo "  sanity:"
head -n2 "$TPM_OUT" | cat -A
awk -F'\t' 'NR==1{print "Header cols:", NF; next} NR==2{print "Data cols:", NF; exit}' "$TPM_OUT"

TPM="$TPM_OUT"

# =========================================
# 5) SUPPA: generate local AS events
# =========================================
echo "[5/8] SUPPA generateEvents..."
python3 "$SUPPA" generateEvents -i "$MERGED" -o "$AS_EVT_DIR/events" -f ioe -e SE SS MX RI FL -p >"$LOGDIR/suppa_generateEvents.log" 2>&1

# =========================================
# 6) SUPPA: PSI per event type
# =========================================
echo "[6/8] SUPPA psiPerEvent..."
for E in "${EVENTS[@]}"; do
  IOE="$AS_EVT_DIR/events_${E}_strict.ioe"
  [ -f "$IOE" ] || { echo "Missing IOE: $IOE"; exit 1; }
  echo "  PSI: $E"
  python3 "$SUPPA" psiPerEvent -i "$IOE" -e "$TPM" -o "$PSI_DIR/${E}" >"$LOGDIR/suppa_psi_${E}.log" 2>&1
done

# =========================================
# helper: compute column indices for sample names (from file header)
# =========================================
get_indices() {
  local file="$1"; shift
  awk -v names="$*" 'BEGIN{FS=OFS="\t";
      n=split(names,a," "); for(i=1;i<=n;i++) want[a[i]]=i}
    NR==1{
      for(i=1;i<=NF;i++) if($i in want) idx[want[$i]]=i
      out="";
      for(i=1;i<=n;i++){
        if(idx[i]) out=out (out?",":"") idx[i];
        else { print "ERROR: sample not in header:", a[i] > "/dev/stderr"; exit 2 }
      }
      print out; exit
    }' "$file"
}

# =========================================
# 7) Split PSI/TPM by condition (by sample names)
# =========================================
echo "[7/8] Split PSI/TPM by condition..."

FR_IDX=$(get_indices "$TPM" "${FR_SAMPLES[@]}")
W_IDX=$(get_indices  "$TPM" "${W_SAMPLES[@]}")

# Split TPM (per-comparison directory)
cut -f1,${FR_IDX} "$TPM" > "$DIFF_DIR/FR.tpm"
cut -f1,${W_IDX}  "$TPM" > "$DIFF_DIR/W.tpm"

# Split PSI per event type (per-comparison directory)
for E in "${EVENTS[@]}"; do
  PSI="$PSI_DIR/${E}.psi"
  [ -f "$PSI" ] || { echo "Missing PSI: $PSI"; exit 1; }
  cut -f1,${FR_IDX} "$PSI" > "$DIFF_DIR/${E}_FR.psi"
  cut -f1,${W_IDX}  "$PSI" > "$DIFF_DIR/${E}_W.psi"
done

# =========================================
# 8) SUPPA: diffSplice (empirical + GC correction)
# =========================================
echo "[8/8] SUPPA diffSplice..."
for E in "${EVENTS[@]}"; do
  IOE="$AS_EVT_DIR/events_${E}_strict.ioe"
  echo "  diffSplice: $E"
  # Output prefix uses comparison name: ${E}_${COMP_NAME}
  python3 "$SUPPA" diffSplice \
    --method empirical -gc \
    --input "$IOE" \
    --psi "$DIFF_DIR/${E}_FR.psi" "$DIFF_DIR/${E}_W.psi" \
    --tpm "$DIFF_DIR/FR.tpm" "$DIFF_DIR/W.tpm" \
    --output "$DIFF_DIR/${E}_${COMP_NAME}" \
    >"$LOGDIR/suppa_diff_${E}_${COMP_NAME}.log" 2>&1
done

echo "✔ All analyses complete.
Comparison: ${COMP_NAME}
- TPM: $TPM
- Events: $AS_EVT_DIR/
- PSI: $PSI_DIR/*.psi
- DiffSplice: $DIFF_DIR/*_${COMP_NAME}.*
"
