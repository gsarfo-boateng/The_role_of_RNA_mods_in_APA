#!/bin/bash

BAM_FILES=(
    "/mnt/george_drive/george/direct_RNAseq/W1118_rep1/W1118_rep1/20250512_1410_MN29498_FBB35082_968b5497/pod5_skip/W1118_rep1_ref.sorted.bam"
    "/mnt/george_drive/george/direct_RNAseq/W1118_rep2/W1118_rep2/20250513_1601_MN29498_FBB35082_95d60745/pod5_skip/W1118_rep2_ref.sorted.bam"
    "/mnt/george_drive/george/direct_RNAseq/FR113N_rep2/FR113N_rep2/FR113N_rep2/20250507_1202_MN29498_FBB33215_0fe8b5ba/pod5_skip/FR113N_05_11_2025_rep2_ref.sorted.bam"
    "/mnt/george_drive/george/direct_RNAseq/FR113N_rep3/FR113N_rep3/20250508_1211_MN29498_FBB33215_84ccd776/pod5_skip/FR113N_05_11_2025_rep3_ref.sorted.bam"
    "/mnt/george_drive/george/direct_RNAseq/W1118_Rep2_mRNA/W1118_Rep2_mRNA/20250703_1210_MN29498_FBB35011_e63c5c87/pod5_skip/W1118_Rep2_mRNA_m6A_m5C_pseu.aligned.sort.bam"
    "/mnt/george_drive/george/direct_RNAseq/W1118_Rep1_mRNA/W1118_Rep1_mRNA/20250702_1515_MN29498_FBB34978_eef47124/pod5_skip/W1118_Rep1_mRNA_m6A_m5C_pseu.aligned.sorted.bam"
    "/mnt/george_drive/george/direct_RNAseq/FR113N_INDURO_mRNA_09_23/FR113N_INDURO_mRNA_09_23/20250923_1413_MN29498_FBC59776_6c8ec548/pod5/FR113N_INDURO_mRNA_09_23.sorted.bam"
    "/mnt/george_drive/george/direct_RNAseq/Mettl3_KO_mRNA_004_rep2/Mettl3_KO_mRNA_004_rep2/20250912_1021_MN29498_FBC38312_0de9666a/pod5/Mettl3_KO_mRNA_rep2.sorted.bam"
    "/mnt/george_drive/george/direct_RNAseq/Mettl3_KO_mRNA_RNA004/Mettl3_KO_mRNA_RNA004/20250910_1505_MN29498_FBC38311_72234379/pod5/Mettl3_KO_mRNA_RNA004.sorted.bam"
    "/mnt/george_drive/george/direct_RNAseq/W1118_Rep3_mRNA/W1118_Rep3_mRNA/20250704_1212_MN29498_FBB35011_b0ca4892/pod5_skip/W1118_Rep3_mRNA_m6A_m5C_pseu.sorted.bam"
    "/mnt/george_drive/george/direct_RNAseq/W1118_induro_mRNA_2025_102/W1118_induro_mRNA_2025_102/20251002_1328_MN29498_FBE06035_a9b4fc34/pod5/W1118_induro_mRNA_2025_102.sorted.bam"
    "/mnt/george_drive/george/direct_RNAseq/w1118_ind_rev_rep1_mRNA/w1118_ind_rev_rep1_mRNA/20250821_1257_MN29498_FBC55934_ad449135/pod5/w1118_ind_rev_rep1_mRNA_m6A_m5C_psi.sorted.bam"
    "/mnt/george_drive/george/direct_RNAseq/w1118_ind_rev_rep1_mRNA_restart/w1118_ind_rev_rep1_mRNA_restart/20250822_1025_MN29498_FBC55934_584462b5/pod5/w1118_ind_rev_rep1_mRNA_restart.sorted.bam"
    "/mnt/george_drive/george/direct_RNAseq/FR113N_Oct_4_mRNA_induro/FR113N_Oct_4_mRNA_induro/20251004_1132_MN29498_FBE00013_f3902e10/pod5/FR113N_Oct_4_mRNA_induro.sorted.bam"
)

OUTPUT_DIR="./modkit_extra_calls"
mkdir -p "$OUTPUT_DIR"

for BAM in "${BAM_FILES[@]}"; do

    SAMPLE_NAME=$(basename "$BAM" .sorted.bam)

    OUT_TSV="${OUTPUT_DIR}/${SAMPLE_NAME}.modkit.tsv"

    echo "Extracting from $BAM --> $OUT_TSV"
#--num-reads 10
    modkit extract calls "$BAM" "$OUT_TSV" --mapped
done
