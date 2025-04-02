#!/bin/bash

A_values=(100 250 500 750 1000 1500) # mMed values
#A_values=(1000 1500)
B_values=(10 20) # mDark values
C_values=(1 100 1000 1500 2000) # ctau values
#C_values=(1500 2000)

for A in "${A_values[@]}"; do
    echo "mMed = ${A}"
    for B in "${B_values[@]}"; do
        echo "mDark = ${B}"
        for C in "${C_values[@]}"; do
	    echo "ctau = ${C}"
	    python3 writeHistograms.py \
		-i filePaths.json \
		-d EMJ_s-channel_mMed-${A}_mDark-${B}_ctau-${C}_unflavored-down \
		-o histograms/histograms_EMJ_s-channel_mMed-${A}_mDark-${B}_ctau-${C}_unflavored-down.root
        done
    done
done
