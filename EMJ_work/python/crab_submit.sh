#!/bin/bash

# Define sample names and corresponding datasets
declare -A samples
samples=(
#    [QCD_Bin-Pt-600to800_Run3Winter25Reco]="/QCD_Bin-Pt-600to800_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
#    [QCD_Bin-Pt-470to600_Run3Winter25Reco]="/QCD_Bin-Pt-470to600_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
#    [QCD_Bin-Pt-30to50_Run3Winter25Reco]="/QCD_Bin-Pt-30to50_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v3/AODSIM"
#    [QCD_Bin-Pt-120to170_Run3Winter25Reco]="/QCD_Bin-Pt-120to170_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
#    [QCD_Bin-Pt-50to80_Run3Winter25Reco]="/QCD_Bin-Pt-50to80_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
#    [QCD_Bin-Pt-170to300_Run3Winter25Reco]="/QCD_Bin-Pt-170to300_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
#    [QCD_Bin-Pt-300to470_Run3Winter25Reco]="/QCD_Bin-Pt-300to470_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
#    [QCD_Bin-Pt-80to120_Run3Winter25Reco]="/QCD_Bin-Pt-80to120_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
#    [QCD_Bin-Pt-800to1000]="/QCD_Bin-Pt-800to1000_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
#    [QCD_Bin-Pt-1000]="/QCD_Bin-Pt-1000_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
    [QCD_PT-15to30]="/QCD_PT-15to30_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM"
    [QCD_PT-30to50]="/QCD_PT-30to50_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM"
#    [QCD_PT-50to80]="/QCD_PT-50to80_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM"
#    [QCD_PT-80to120]="/QCD_PT-80to120_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM"
#    [QCD_PT-120to170]="/QCD_PT-120to170_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM"
#    [QCD_PT-170to300]="/QCD_PT-170to300_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM"
#    [QCD_PT-300to470]="/QCD_PT-300to470_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM"
#    [QCD_PT-470to600]="/QCD_PT-470to600_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM"
#    [QCD_PT-600to800]="/QCD_PT-600to800_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM"
#    [QCD_PT-800to1000]="/QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM"
#    [QCD_PT-1000to1400]="/QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM"
#    [QCD_PT-1400to1800]="/QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM"
#    [QCD_PT-1800to2400]="/QCD_PT-1800to2400_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM"
#    [QCD_PT-2400to3200]="/QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM"
#    [QCD_PT-3200]="/QCD_PT-3200_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM"
)

todaysDate=$(date +%Y%m%d)

for sample_name in "${!samples[@]}"; do
    sample_dataset="${samples[$sample_name]}"

    # Create a CRAB config file
    cat > crab_cfg_$sample_name.py <<EOL
import os
import datetime
from CRABClient.UserUtilities import config
config = config()

todaysDate = "$todaysDate"

config.General.requestName = f'Run3EMJJetStudy_${sample_name}_${todaysDate}'
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = '/afs/cern.ch/user/e/ekauffma/crabWorkArea'

config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = f'${CMSSW_BASE}/src/EMJ/EMJ_work/python/makeNtuples_cfg.py'
config.JobType.pluginName = 'Analysis'
config.Data.outLFNDirBase = '/store/group/lpcemj/EMJAnalysis2025/'

config.JobType.maxMemoryMB = 5000

config.Data.inputDataset = "$sample_dataset"
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 10000
config.Data.publication = False
config.Data.outputDatasetTag = f'Run3EMJJetStudy_${sample_name}_${todaysDate}'
config.Data.partialDataset = True

config.Site.storageSite = 'T3_US_FNALLPC'
EOL

    # Submit the job
    echo "Submitting job for $sample_name with dataset $sample_dataset"
    crab submit -c crab_cfg_$sample_name.py

done

