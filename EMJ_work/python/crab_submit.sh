#!/bin/bash

# Define sample names and corresponding datasets
declare -A samples
samples=(
    [QCD_Bin-Pt-600to800_Run3Winter25Reco]="/QCD_Bin-Pt-600to800_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
    [QCD_Bin-Pt-470to600_Run3Winter25Reco]="/QCD_Bin-Pt-470to600_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
    [QCD_Bin-Pt-30to50_Run3Winter25Reco]="/QCD_Bin-Pt-30to50_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v3/AODSIM"
    [QCD_Bin-Pt-120to170_Run3Winter25Reco]="/QCD_Bin-Pt-120to170_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
    [QCD_Bin-Pt-50to80_Run3Winter25Reco]="/QCD_Bin-Pt-50to80_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
    [QCD_Bin-Pt-170to300_Run3Winter25Reco]="/QCD_Bin-Pt-170to300_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
    [QCD_Bin-Pt-300to470_Run3Winter25Reco]="/QCD_Bin-Pt-300to470_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
    [QCD_Bin-Pt-80to120_Run3Winter25Reco]="/QCD_Bin-Pt-80to120_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
    [QCD_Bin-Pt-800to1000]="/QCD_Bin-Pt-800to1000_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
    [QCD_Bin-Pt-1000]="/QCD_Bin-Pt-1000_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-142X_mcRun3_2025_realistic_v7-v2/AODSIM"
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
config.Data.outLFNDirBase = '/store/user/ekauffma/'

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

