#!/bin/bash

# Loop through all BAM files in the current directory
for bam in *.bam; do
  # Split the current BAM file into plus and minus strands
  samtools view -F 16 -b -o ${bam%.bam}_plus.bam $bam 
  samtools view -f 16 -b -o ${bam%.bam}_minus.bam $bam
  #samtools view -h $bam | awk '$2 ~ /^[-+]$/ {print $0}' | samtools view -b - > ${bam%.bam}_minus.bam
  
  # Calculate the total number of reads in the plus strand BAM file
  total_reads_plus=$(samtools view -c ${bam%.bam}_plus.bam)
  # Divide the total number of reads by 1000000
  scale_plus=$(echo "scale=6; $total_reads_plus/1000000" | bc)
  # Divide the plus strand BAM file by the total number of reads and convert to BedGraph
  bedtools genomecov -ibam ${bam%.bam}_plus.bam -bg -scale $scale_plus > ${bam%.bam}_plus.bedgraph
  
  # Calculate the total number of reads in the minus strand BAM file
  total_reads_minus=$(samtools view -c ${bam%.bam}_minus.bam)
  # Divide the total number of reads by 1000000
  scale_minus=$(echo "scale=6; $total_reads_minus/1000000" | bc)
  # Divide the minus strand BAM file by the total number of reads and convert to BedGraph
  bedtools genomecov -ibam ${bam%.bam}_minus.bam -bg -scale $scale_minus > ${bam%.bam}_minus.bedgraph
  
  #sort the resulting bedgraphs
  sort -k1,1 -k2,2n ${bam%.bam}_plus.bedgraph > ${bam%.bam}_plus_sorted.bedgraph
  sort -k1,1 -k2,2n ${bam%.bam}_minus.bedgraph > ${bam%.bam}_minus_sorted.bedgraph
  
  #convert the minus bedgraphs with awk '{ $4 *= -1 } 1'
  awk '{ $4 *= -1 } 1' ${bam%.bam}_minus_sorted.bedgraph > ${bam%.bam}_minus_sorted_converted.bedgraph
  
  # Convert the resulting BedGraph files into BigWig files
  bedGraphToBigWig ${bam%.bam}_plus_sorted.bedgraph /media/judy/george/m6a/dm6.chrom.sizes ${bam%.bam}_plus.bw
  bedGraphToBigWig ${bam%.bam}_minus_sorted_converted.bedgraph /media/judy/george/m6a//dm6.chrom.sizes ${bam%.bam}_minus.bw
done
