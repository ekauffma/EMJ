#!/bin/bash

A_values=(100 250 500 750 1000 1500)
B_values=(10 20)
C_values=(1 100 1000 1500 2000)

for A in "${A_values[@]}"; do
    for B in "${B_values[@]}"; do
         for C in "${C_values[@]}"; do
             echo "mMed  = ${A}, mDark = ${B}, ctau = ${C}"
	     python3 plot2DHistograms.py \
                 -i "histograms/histograms_EMJ_s-channel_mMed-${A}_mDark-${B}_ctau-${C}_unflavored-down.root" \
                 -d "EMJ_s-channel_mMed-${A}_mDark-${B}_ctau-${C}_unflavored-down" \
                 -o plots/	
	 done	
    done
done
