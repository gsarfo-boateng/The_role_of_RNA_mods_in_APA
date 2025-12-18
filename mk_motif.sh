#!/usr/bin/env bash
ref="/mnt/george_drive/george/m6a/dm6.fasta"
mkdir -p motifs logs

for i in *.bed; do
  base="${i%.bed}"
  modkit find-motifs -i "$i" -r "$ref" -o "motifs/${base}.motifs.tsv" --threads 32 --log "logs/${base}.log"
done
