#!/bin/bash

QCD_Bins=("15to30" "30to50" "50to80" "80to120" "120to170" "170to300" "300to470" "470to600" "600to800" "800to1000" "1000to1400" "1400to1800" "1800to2400" "2400to3200" "3200") # QCD bins

for Bin in "${QCD_Bins[@]}"; do
    echo "PT Range = ${Bin}"
    python3 writeHistograms.py \
	-i filePaths.json \
	-d QCD_PT-${Bin} \
	-o histograms/histograms_QCD_Bin-Pt-${Bin}.root
    done
