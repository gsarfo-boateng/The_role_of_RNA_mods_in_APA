#!/usr/bin/env bash

REF="/mnt/raid1/george/direct_RNA/"
THREADS=8

# Replicates
# FR113N (3 reps)
FR113N_FILES=(
  "FR113N_05_11_2025_rep2_ref.sorted.bed.gz"
  "FR113N_05_11_2025_rep3_ref.sorted.bed.gz"
  "FR113N_Oct_4_mRNA_induro.sorted.bed.gz")

# W1118 (3 reps)
W1118_FILES=(
  "W1118_induro_mRNA_2025_102.sorted.bed.gz"
  "W1118_Rep1_mRNA_m6A_m5C_pseu.aligned.sorted.bed.gz"
  "W1118_rep2_ref.sorted.bed.gz")

# Mettl3-KO (2 reps)
KO_FILES=(
  "Mettl3_KO_mRNA_RNA004.sorted.bed.gz"
  "Mettl3_KO_mRNA_rep2.sorted.bed.gz")

# Helpers
check_files() {
  for f in "$@"; do
    [[ -e "$f" ]] || { echo "[ERROR] Missing: $f" >&2; exit 1; }
  done
}

run_dmr_to_file() {
  local base="$1"
  local outfile="$2"
  shift 2
  local -a a_files=("$1"); shift
  local -a b_files=("$1"); shift

  echo "[INFO] base=${base} -> ${outfile}"
  modkit dmr pair $(printf ' -a %q' "${a_files[@]}")  $(printf ' -b %q' "${b_files[@]}") --ref "$REF" --base "$base" --out-path "$outfile" -t "$THREADS" -f --header
}

# sanity checks
check_files "$REF"
check_files "${FR113N_FILES[@]}" "${W1118_FILES[@]}" "${KO_FILES[@]}"

# CONTRAST 1: FR113N vs Mettl3-KO
# Outputs: FR113NOct4Mettl3KOdmr_{a,C,T}.bed
run_dmr_to_file A "FR113NOct4Mettl3KOdmr_a.bed" "${FR113N_FILES[@]}" "${KO_FILES[@]}"
run_dmr_to_file C "FR113NOct4Mettl3KOdmr_C.bed" "${FR113N_FILES[@]}" "${KO_FILES[@]}"
run_dmr_to_file T "FR113NOct4Mettl3KOdmr_T.bed" "${FR113N_FILES[@]}" "${KO_FILES[@]}"

# CONTRAST 2: W1118 vs Mettl3-KO
# Outputs: w1118_vs_Mettl3KOdmr_{a,C,T}.bed
run_dmr_to_file A "w1118_vs_Mettl3KOdmr_a.bed" "${W1118_FILES[@]}" "${KO_FILES[@]}"
run_dmr_to_file C "w1118_vs_Mettl3KOdmr_C.bed" "${W1118_FILES[@]}" "${KO_FILES[@]}"
run_dmr_to_file T "w1118_vs_Mettl3KOdmr_T.bed" "${W1118_FILES[@]}" "${KO_FILES[@]}"

# CONTRAST 3: FR113N vs W1118
# Outputs: FR113NOct4_05_11_W1118_Rep2_mRNA_m6_W1118_induro_mRNA_2025_102_dmr_{a,C,T}.bed
run_dmr_to_file A "FR113NOct4_05_11_W1118_Rep2_mRNA_m6_W1118_induro_mRNA_2025_102_dmr_a.bed" "${FR113N_FILES[@]}" "${W1118_FILES[@]}"
run_dmr_to_file C "FR113NOct4_05_11_W1118_Rep2_mRNA_m6_W1118_induro_mRNA_2025_102_dmr_C.bed" "${FR113N_FILES[@]}" "${W1118_FILES[@]}"
run_dmr_to_file T "FR113NOct4_05_11_W1118_Rep2_mRNA_m6_W1118_induro_mRNA_2025_102_dmr_T.bed" "${FR113N_FILES[@]}" "${W1118_FILES[@]}"

echo "[DONE] DMR .bed runs complete."
