#!/bin/bash

A_values=(100 250 500 750 1000 1500) # mMed values
B_values=(10 20) # mDark values
C_values=(1 100 1000) # ctau values

for A in "${A_values[@]}"; do
    for B in "${B_values[@]}"; do
        for C in "${C_values[@]}"; do
	    python3 writeHistograms.py \
		-i filePaths.json \
		-d EMJ_s-channel_mMed-${A}_mDark-${B}_ctau-${C}_unflavored-down \
		-o histograms/histograms_EMJ_s-channel_mMed-${A}_mDark-${B}_ctau-${C}_unflavored-down.root
        done
    done
done
