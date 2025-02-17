#!/bin/bash
/bin/pwd
/bin/uname -a
/bin/hostname

source /cvmfs/cms.cern.ch/cmsset_default.sh
export X509_USER_PROXY=$PWD/x509_proxy

voms-proxy-info --all

cmsrel CMSSW_14_2_2

cd CMSSW_14_2_2/src
mkdir EMJ
cd EMJ

cp -r ../../../EMJ_work .

cmsenv
scram b -j 8

echo "finished setting up cmssw"

cd EMJ_work/python

cmsRun makeNtuples_cfg.py inputFiles=$1

outputFileName=$(basename "$1" | sed 's/^step_RECO_//')
xrdcp output.root root://cmseos.fnal.gov//store/user/ekauffma/EMJ/signal_ntuples/ntuple_$outputFileName
rm output.root
