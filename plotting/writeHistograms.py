###################################################################################################
#   writeHistograms.py                                                                            #
#   Description: create histograms of jet constituent info                                        #
#   Author: Elliott Kauffman                                                                      #
###################################################################################################

###################################################################################################
# IMPORTS

import argparse
#from distributed import Client
import json
import multiprocessing
import ROOT
from time import time

###################################################################################################
# RUN OPTIONS

treesForDF = [
    "JetInfoNtuplizer/PFCand",
]

###################################################################################################
# METHODS

def parseArgs() -> argparse.Namespace:

    p = argparse.ArgumentParser()
    p.add_argument(
        "--input",
        "-i",
        help="Name of the json file containing the filepaths to process.",
        default="filePathsTest.json",
    )
    p.add_argument(
        "--dataset",
        "-d",
        help="Name of json key for dataset we want to process",
    )
    p.add_argument(
        "--output",
        "-o",
        help="Name of the file where histograms will be stored. If it already exists, contents are overwritten.",
        default="histograms.root",
    )
    return p.parse_args()


# createChain: returns ROOT TChain of the input fileList, using all trees in treeList
def createChain(
    fileList,
    treeList
) -> ROOT.TChain :

    chainList = []
    for treeName in treeList:
        chainList.append(ROOT.TChain(treeName))
        for filePath in fileList:
            chainList[-1].Add(filePath)

    returnChain = chainList[0]
    if len(chainList)>1:
        for chain in chainList[1:]:
            returnChain.AddFriend(chain)

    return returnChain

def declareCppMethods():

    # function for applying jet mask to nested variable
    ROOT.gInterpreter.Declare(R"""
#include <vector>
#include "ROOT/RVec.hxx"

using namespace ROOT::VecOps;

// Function to filter jets based on a mask
RVec<RVec<double>> filterJetsByMask(const RVec<RVec<double>>& jetData, const RVec<int>& mask) {
    RVec<RVec<double>> filteredJets;
    for (size_t i = 0; i < jetData.size(); ++i) {
        if (mask[i] == 1) {  // Keep only jets that pass the mask condition
            filteredJets.push_back(jetData[i]);
        }
    }
    return filteredJets;
}
""")

    # function for computing averages of constituent info per jet
    ROOT.gInterpreter.Declare("""
#include <vector>
#include "ROOT/RVec.hxx"

RVec<double> computeJetAverages(const RVec<RVec<double>>& perJetData) {
    RVec<double> jetAverages;
    for (const auto& jetConstituents : perJetData) {
        if (!jetConstituents.empty()) {
            double sum = 0;
            for (double value : jetConstituents) {
                sum+=value;
            }
            jetAverages.push_back(sum/jetConstituents.size());
        }
        else jetAverages.push_back(0.0);
    }
    return jetAverages;
}
""")

    # function for computing pT of jet constituents
    ROOT.gInterpreter.Declare("""
#include <vector>
#include "ROOT/RVec.hxx"

RVec<RVec<double>> computeJetConstituentPt(
    const RVec<RVec<double>>& JetConstituentPx,
    const RVec<RVec<double>>& JetConstituentPy)
{
    RVec<RVec<double>> JetConstituentPt;
    for (size_t i = 0; i < JetConstituentPx.size(); ++i) {
        RVec<double> ptValues;
        for (size_t j = 0; j < JetConstituentPx[i].size(); ++j) {
            double pt = std::sqrt(JetConstituentPx[i][j] * JetConstituentPx[i][j] + JetConstituentPy[i][j] * JetConstituentPy[i][j]);
            ptValues.push_back(pt);
        }
        JetConstituentPt.push_back(ptValues);
    }
    return JetConstituentPt;
}
""")

    # function for computing pt-weighted averages of constituent info per jet
    ROOT.gInterpreter.Declare("""
#include <vector>
#include "ROOT/RVec.hxx"

RVec<double> computePtWeightedJetAverages(
    const RVec<RVec<double>>& perJetConstituentVar,
    const RVec<RVec<double>>& perJetConstituentPt)
{
        RVec<double> jetWeightedAverages;
    for(size_t i = 0; i < perJetConstituentVar.size(); ++i) {
        double sum = 0.0;
        double ptsum = 0.0;
        for (size_t j = 0; j < perJetConstituentVar[i].size(); ++j) {
            sum += perJetConstituentPt[i][j] * perJetConstituentVar[i][j];
            ptsum += perJetConstituentPt[i][j];
        }
        if (ptsum!=0.0) jetWeightedAverages.push_back(sum / ptsum);
        else jetWeightedAverages.push_back(0.0);
    }
    
    return jetWeightedAverages;
}
""")

    # function for computing medians of constituent info per jet
    ROOT.gInterpreter.Declare("""
#include <vector>
#include <algorithm>
#include "ROOT/RVec.hxx"
    
RVec<double> computeJetMedians(const RVec<RVec<double>>& perJetConstituentVar) {
    RVec<double> jetMedians;
    for(size_t i = 0; i < perJetConstituentVar.size(); ++i) {
        if (!perJetConstituentVar[i].empty()) {
            RVec<double> sortedData = perJetConstituentVar[i];
            std::sort(sortedData.begin(), sortedData.end());
            size_t n = sortedData.size();
            if (n % 2 == 0) {
                jetMedians.push_back((sortedData[n / 2 - 1] + sortedData[n / 2]) / 2.0);
            }
            else jetMedians.push_back(sortedData[n / 2]);
        }
        else jetMedians.push_back(0.0);
    }
    return jetMedians;
}
""")

    # function for computing minimum of constituent info per jet
    ROOT.gInterpreter.Declare("""
#include <vector>
#include <algorithm>
#include "ROOT/RVec.hxx"
    
RVec<double> computeJetMinimums(const RVec<RVec<double>>& perJetConstituentVar) {
    RVec<double> jetMinimums;
    for(size_t i = 0; i < perJetConstituentVar.size(); ++i) {
        double minimum = 1.0;
        for (size_t j = 0; j < perJetConstituentVar[i].size(); ++j) {
            if (perJetConstituentVar[i][j] < minimum) minimum = perJetConstituentVar[i][j];
        }
        jetMinimums.push_back(minimum);
    }
    return jetMinimums;
}
""")

    # function for computing maximum of constituent info per jet
    ROOT.gInterpreter.Declare("""
#include <vector>
#include <algorithm>
#include "ROOT/RVec.hxx"
    
RVec<double> computeJetMaximums(const RVec<RVec<double>>& perJetConstituentVar) {
    RVec<double> jetMaximums;
    for(size_t i = 0; i < perJetConstituentVar.size(); ++i) {
        double maximum = 0.0;
        for (size_t j = 0; j < perJetConstituentVar[i].size(); ++j) {
            if (perJetConstituentVar[i][j] > maximum) maximum = perJetConstituentVar[i][j];
        }
        jetMaximums.push_back(maximum);
    }
    return jetMaximums;
}
""")
    
    # function for filtering out constituents that leave no energy in the HCAL
    ROOT.gInterpreter.Declare("""
#include <vector>
#include <cmath>
using namespace ROOT::VecOps;

RVec<int> getValidConstituentIndices(
    const RVec<float>& hcalDepthEF1,
    const RVec<float>& hcalDepthEF2,
    const RVec<float>& hcalDepthEF3,
    const RVec<float>& hcalDepthEF4,
    const RVec<float>& hcalDepthEF5,
    const RVec<float>& hcalDepthEF6,
    const RVec<float>& hcalDepthEF7) {

    RVec<int> valid_indices;
    for (size_t i = 0; i < hcalDepthEF1.size(); ++i) {
        bool all_zero = (hcalDepthEF1[i] == 0) && (hcalDepthEF2[i] == 0) &&
                        (hcalDepthEF3[i] == 0) && (hcalDepthEF4[i] == 0) &&
                        (hcalDepthEF5[i] == 0) && (hcalDepthEF6[i] == 0) &&
                        (hcalDepthEF7[i] == 0);
        if (!all_zero) valid_indices.push_back(i);
    }
    return valid_indices;
}
""")

    # function to filter a vector based on valid indices
    ROOT.gInterpreter.Declare("""
#include <vector>
#include <cmath>
using namespace ROOT::VecOps;

template <typename T>
RVec<T> filterVector(const RVec<T>& vec, const RVec<int>& indices) {
    RVec<T> filtered;
    for (int idx : indices) {
        filtered.push_back(vec[idx]);
    }
    return filtered;
}
    
""")

    # function to recompute jet number of constituents after filtering
    ROOT.gInterpreter.Declare("""
#include <vector>
#include <cmath>
using namespace ROOT::VecOps;

RVec<int> recomputeJetNConstituents(const RVec<int>& original_nConstituents, const RVec<int>& valid_indices) {
    RVec<int> new_nConstituents;
    int index = 0;
    for (auto n : original_nConstituents) {
        int count = 0;
        for (int i = 0; i < n; i++) {
            if (index < valid_indices.size() && valid_indices[index] == i) {
                count++;
                index++;
            }
        }
        new_nConstituents.push_back(count);
    }
    return new_nConstituents;
}
    
""")
    


    # function to reconstruct per-jet structure after filtering
    ROOT.gInterpreter.Declare("""
#include <vector>
#include <cmath>
using namespace ROOT::VecOps;
    
template <typename T>
RVec<RVec<T>> reconstructPerJet(const RVec<T>& vec, const RVec<int>& Jet_nConstituent) {
    RVec<RVec<T>> perJet;
    int offset = 0;
    for (size_t i = 0; i < Jet_nConstituent.size(); ++i) {
        int n = Jet_nConstituent[i];
        perJet.emplace_back(vec.begin() + offset, vec.begin() + offset + n);
        offset += n;
    }
    return perJet;
}
    
""")

    return
    

def prepareDataframe(df):

    # include events with at least 1 jet with pT > 20 GeV
    df = df.Define("Jet_mask", "(Jet_Pt > 20) && (abs(Jet_Eta)<2.5)") # apply jet selections here
    df = df.Filter("Sum(Jet_mask)>=1") # require at least one good jet per event
    print("Filtered events to contain at least one good jet")
    
    # filter constituents
    df = df.Define("valid_indices",
                   "getValidConstituentIndices(Constituent_hcalDepthEF1, Constituent_hcalDepthEF2, Constituent_hcalDepthEF3, Constituent_hcalDepthEF4, Constituent_hcalDepthEF5, Constituent_hcalDepthEF6, Constituent_hcalDepthEF7)"
                   )
                   
    df = df.Define("Jet_nFilteredConstituent", "recomputeJetNConstituents(Jet_nConstituent, valid_indices)")
                   
    df = df.Define("filteredConstituent_Px", "filterVector(Constituent_Px, valid_indices)")
    df = df.Define("filteredConstituent_Py", "filterVector(Constituent_Py, valid_indices)")
    df = df.Define("filteredConstituent_Pz", "filterVector(Constituent_Pz, valid_indices)")
    df = df.Define("filteredConstituent_E", "filterVector(Constituent_E, valid_indices)")
    
    df = df.Define("filteredConstituent_hcalDepthEF1", "filterVector(Constituent_hcalDepthEF1, valid_indices)")
    df = df.Define("filteredConstituent_hcalDepthEF2", "filterVector(Constituent_hcalDepthEF2, valid_indices)")
    df = df.Define("filteredConstituent_hcalDepthEF3", "filterVector(Constituent_hcalDepthEF3, valid_indices)")
    df = df.Define("filteredConstituent_hcalDepthEF4", "filterVector(Constituent_hcalDepthEF4, valid_indices)")
    df = df.Define("filteredConstituent_hcalDepthEF5", "filterVector(Constituent_hcalDepthEF5, valid_indices)")
    df = df.Define("filteredConstituent_hcalDepthEF6", "filterVector(Constituent_hcalDepthEF6, valid_indices)")
    df = df.Define("filteredConstituent_hcalDepthEF7", "filterVector(Constituent_hcalDepthEF7, valid_indices)")
    print("Filtered out constituents with no deposits in HCAL")
    
    
    # unflatten jet constituent info
    df = df.Define("Jet_Constituent_Px", "reconstructPerJet(filteredConstituent_Px, Jet_nFilteredConstituent)")
    df = df.Define("Jet_Constituent_Py", "reconstructPerJet(filteredConstituent_Py, Jet_nFilteredConstituent)")
    df = df.Define("Jet_Constituent_Pz", "reconstructPerJet(filteredConstituent_Pz, Jet_nFilteredConstituent)")
    df = df.Define("Jet_Constituent_E", "reconstructPerJet(filteredConstituent_E, Jet_nFilteredConstituent)")
    print("Unflatted jet constituent kinematics variables")
    
    df = df.Define("Jet_Constituent_hcalDepthEF1", "reconstructPerJet(filteredConstituent_hcalDepthEF1, Jet_nFilteredConstituent)")
    df = df.Define("Jet_Constituent_hcalDepthEF2", "reconstructPerJet(filteredConstituent_hcalDepthEF2, Jet_nFilteredConstituent)")
    df = df.Define("Jet_Constituent_hcalDepthEF3", "reconstructPerJet(filteredConstituent_hcalDepthEF3, Jet_nFilteredConstituent)")
    df = df.Define("Jet_Constituent_hcalDepthEF4", "reconstructPerJet(filteredConstituent_hcalDepthEF4, Jet_nFilteredConstituent)")
    df = df.Define("Jet_Constituent_hcalDepthEF5", "reconstructPerJet(filteredConstituent_hcalDepthEF5, Jet_nFilteredConstituent)")
    df = df.Define("Jet_Constituent_hcalDepthEF6", "reconstructPerJet(filteredConstituent_hcalDepthEF6, Jet_nFilteredConstituent)")
    df = df.Define("Jet_Constituent_hcalDepthEF7", "reconstructPerJet(filteredConstituent_hcalDepthEF7, Jet_nFilteredConstituent)")
    print("Unflatted jet constituent HCAL depth fractions")
    
    # filter out bad jets
    df = df.Define("goodJet_Pt", "Jet_Pt[Jet_mask]")
    df = df.Define("goodJet_Eta", "Jet_Eta[Jet_mask]")
    df = df.Define("goodJet_Phi", "Jet_Phi[Jet_mask]")
    df = df.Define("goodJet_E", "Jet_E[Jet_mask]")
    df = df.Define("goodJet_nFilteredConstituent", "Jet_nFilteredConstituent[Jet_mask]")
    df = df.Define("goodJet_Constituent_Px", "filterJetsByMask(Jet_Constituent_Px, Jet_mask)")
    df = df.Define("goodJet_Constituent_Py", "filterJetsByMask(Jet_Constituent_Py, Jet_mask)")
    df = df.Define("goodJet_Constituent_Pz", "filterJetsByMask(Jet_Constituent_Pz, Jet_mask)")
    df = df.Define("goodJet_Constituent_E", "filterJetsByMask(Jet_Constituent_E, Jet_mask)")
    df = df.Define("goodJet_Constituent_hcalDepthEF1", "filterJetsByMask(Jet_Constituent_hcalDepthEF1, Jet_mask)")
    df = df.Define("goodJet_Constituent_hcalDepthEF2", "filterJetsByMask(Jet_Constituent_hcalDepthEF2, Jet_mask)")
    df = df.Define("goodJet_Constituent_hcalDepthEF3", "filterJetsByMask(Jet_Constituent_hcalDepthEF3, Jet_mask)")
    df = df.Define("goodJet_Constituent_hcalDepthEF4", "filterJetsByMask(Jet_Constituent_hcalDepthEF4, Jet_mask)")
    df = df.Define("goodJet_Constituent_hcalDepthEF5", "filterJetsByMask(Jet_Constituent_hcalDepthEF5, Jet_mask)")
    df = df.Define("goodJet_Constituent_hcalDepthEF6", "filterJetsByMask(Jet_Constituent_hcalDepthEF6, Jet_mask)")
    df = df.Define("goodJet_Constituent_hcalDepthEF7", "filterJetsByMask(Jet_Constituent_hcalDepthEF7, Jet_mask)")
    print("Filtered out bad jets")
    
    # compute averages of jet constituent info per jet
    df = df.Define("goodJet_avgConstituentHcalDepthEF1", "computeJetAverages(goodJet_Constituent_hcalDepthEF1)")
    df = df.Define("goodJet_avgConstituentHcalDepthEF2", "computeJetAverages(goodJet_Constituent_hcalDepthEF2)")
    df = df.Define("goodJet_avgConstituentHcalDepthEF3", "computeJetAverages(goodJet_Constituent_hcalDepthEF3)")
    df = df.Define("goodJet_avgConstituentHcalDepthEF4", "computeJetAverages(goodJet_Constituent_hcalDepthEF4)")
    df = df.Define("goodJet_avgConstituentHcalDepthEF5", "computeJetAverages(goodJet_Constituent_hcalDepthEF5)")
    df = df.Define("goodJet_avgConstituentHcalDepthEF6", "computeJetAverages(goodJet_Constituent_hcalDepthEF6)")
    df = df.Define("goodJet_avgConstituentHcalDepthEF7", "computeJetAverages(goodJet_Constituent_hcalDepthEF7)")
    print("Computed averages of constituent HCAL depth EF per jet")
    
    # compute jet constituent Pt
    df = df.Define("goodJet_Constituent_Pt", "computeJetConstituentPt(goodJet_Constituent_Px, goodJet_Constituent_Py)")
    print("Defined jet constituent transverse momentum")
    
    # compute pt-weighted averages of jet constituent info per jet
    df = df.Define("goodJet_pTWeightedAvgConstituentHcalDepthEF1", "computePtWeightedJetAverages(goodJet_Constituent_hcalDepthEF1, goodJet_Constituent_Pt)")
    df = df.Define("goodJet_pTWeightedAvgConstituentHcalDepthEF2", "computePtWeightedJetAverages(goodJet_Constituent_hcalDepthEF2, goodJet_Constituent_Pt)")
    df = df.Define("goodJet_pTWeightedAvgConstituentHcalDepthEF3", "computePtWeightedJetAverages(goodJet_Constituent_hcalDepthEF3, goodJet_Constituent_Pt)")
    df = df.Define("goodJet_pTWeightedAvgConstituentHcalDepthEF4", "computePtWeightedJetAverages(goodJet_Constituent_hcalDepthEF4, goodJet_Constituent_Pt)")
    df = df.Define("goodJet_pTWeightedAvgConstituentHcalDepthEF5", "computePtWeightedJetAverages(goodJet_Constituent_hcalDepthEF5, goodJet_Constituent_Pt)")
    df = df.Define("goodJet_pTWeightedAvgConstituentHcalDepthEF6", "computePtWeightedJetAverages(goodJet_Constituent_hcalDepthEF6, goodJet_Constituent_Pt)")
    df = df.Define("goodJet_pTWeightedAvgConstituentHcalDepthEF7", "computePtWeightedJetAverages(goodJet_Constituent_hcalDepthEF7, goodJet_Constituent_Pt)")
    print("Computed pt-weighted averages of constituent HCAL depth EF per jet")
    
    # compute medians of jet constituent info per jet
    df = df.Define("goodJet_medConstituentHcalDepthEF1", "computeJetMedians(goodJet_Constituent_hcalDepthEF1)")
    df = df.Define("goodJet_medConstituentHcalDepthEF2", "computeJetMedians(goodJet_Constituent_hcalDepthEF2)")
    df = df.Define("goodJet_medConstituentHcalDepthEF3", "computeJetMedians(goodJet_Constituent_hcalDepthEF3)")
    df = df.Define("goodJet_medConstituentHcalDepthEF4", "computeJetMedians(goodJet_Constituent_hcalDepthEF4)")
    df = df.Define("goodJet_medConstituentHcalDepthEF5", "computeJetMedians(goodJet_Constituent_hcalDepthEF5)")
    df = df.Define("goodJet_medConstituentHcalDepthEF6", "computeJetMedians(goodJet_Constituent_hcalDepthEF6)")
    df = df.Define("goodJet_medConstituentHcalDepthEF7", "computeJetMedians(goodJet_Constituent_hcalDepthEF7)")
    print("Computed medians of constituent HCAL depth EF per jet")
    
    # compute minimums of jet constituent info per jet
    df = df.Define("goodJet_minConstituentHcalDepthEF1", "computeJetMinimums(goodJet_Constituent_hcalDepthEF1)")
    df = df.Define("goodJet_minConstituentHcalDepthEF2", "computeJetMinimums(goodJet_Constituent_hcalDepthEF2)")
    df = df.Define("goodJet_minConstituentHcalDepthEF3", "computeJetMinimums(goodJet_Constituent_hcalDepthEF3)")
    df = df.Define("goodJet_minConstituentHcalDepthEF4", "computeJetMinimums(goodJet_Constituent_hcalDepthEF4)")
    df = df.Define("goodJet_minConstituentHcalDepthEF5", "computeJetMinimums(goodJet_Constituent_hcalDepthEF5)")
    df = df.Define("goodJet_minConstituentHcalDepthEF6", "computeJetMinimums(goodJet_Constituent_hcalDepthEF6)")
    df = df.Define("goodJet_minConstituentHcalDepthEF7", "computeJetMinimums(goodJet_Constituent_hcalDepthEF7)")
    print("Computed minimums of constituent HCAL depth EF per jet")
    
    # compute maximums of jet constituent info per jet
    df = df.Define("goodJet_maxConstituentHcalDepthEF1", "computeJetMaximums(goodJet_Constituent_hcalDepthEF1)")
    df = df.Define("goodJet_maxConstituentHcalDepthEF2", "computeJetMaximums(goodJet_Constituent_hcalDepthEF2)")
    df = df.Define("goodJet_maxConstituentHcalDepthEF3", "computeJetMaximums(goodJet_Constituent_hcalDepthEF3)")
    df = df.Define("goodJet_maxConstituentHcalDepthEF4", "computeJetMaximums(goodJet_Constituent_hcalDepthEF4)")
    df = df.Define("goodJet_maxConstituentHcalDepthEF5", "computeJetMaximums(goodJet_Constituent_hcalDepthEF5)")
    df = df.Define("goodJet_maxConstituentHcalDepthEF6", "computeJetMaximums(goodJet_Constituent_hcalDepthEF6)")
    df = df.Define("goodJet_maxConstituentHcalDepthEF7", "computeJetMaximums(goodJet_Constituent_hcalDepthEF7)")
    print("Computed maximums of constituent HCAL depth EF per jet")
    
    # define leading jet values
    df = df.Define("leadJet_Pt", "goodJet_Pt[0]")
    df = df.Define("leadJet_Eta", "goodJet_Eta[0]")
    df = df.Define("leadJet_Phi", "goodJet_Phi[0]")
    df = df.Define("leadJet_E", "goodJet_E[0]")
    df = df.Define("leadJet_nConstituent", "goodJet_nFilteredConstituent[0]")
    
    # filter out events where the leading jet has no remaining constituents
    df = df.Filter("leadJet_nConstituent > 0")
    
    df = df.Define("leadJet_avgConstituentHcalDepthEF1", "goodJet_avgConstituentHcalDepthEF1[0]")
    df = df.Define("leadJet_avgConstituentHcalDepthEF2", "goodJet_avgConstituentHcalDepthEF2[0]")
    df = df.Define("leadJet_avgConstituentHcalDepthEF3", "goodJet_avgConstituentHcalDepthEF3[0]")
    df = df.Define("leadJet_avgConstituentHcalDepthEF4", "goodJet_avgConstituentHcalDepthEF4[0]")
    df = df.Define("leadJet_avgConstituentHcalDepthEF5", "goodJet_avgConstituentHcalDepthEF5[0]")
    df = df.Define("leadJet_avgConstituentHcalDepthEF6", "goodJet_avgConstituentHcalDepthEF6[0]")
    df = df.Define("leadJet_avgConstituentHcalDepthEF7", "goodJet_avgConstituentHcalDepthEF7[0]")
    
    df = df.Define("leadJet_pTWeightedAvgConstituentHcalDepthEF1", "goodJet_pTWeightedAvgConstituentHcalDepthEF1[0]")
    df = df.Define("leadJet_pTWeightedAvgConstituentHcalDepthEF2", "goodJet_pTWeightedAvgConstituentHcalDepthEF2[0]")
    df = df.Define("leadJet_pTWeightedAvgConstituentHcalDepthEF3", "goodJet_pTWeightedAvgConstituentHcalDepthEF3[0]")
    df = df.Define("leadJet_pTWeightedAvgConstituentHcalDepthEF4", "goodJet_pTWeightedAvgConstituentHcalDepthEF4[0]")
    df = df.Define("leadJet_pTWeightedAvgConstituentHcalDepthEF5", "goodJet_pTWeightedAvgConstituentHcalDepthEF5[0]")
    df = df.Define("leadJet_pTWeightedAvgConstituentHcalDepthEF6", "goodJet_pTWeightedAvgConstituentHcalDepthEF6[0]")
    df = df.Define("leadJet_pTWeightedAvgConstituentHcalDepthEF7", "goodJet_pTWeightedAvgConstituentHcalDepthEF7[0]")
    
    df = df.Define("leadJet_medConstituentHcalDepthEF1", "goodJet_medConstituentHcalDepthEF1[0]")
    df = df.Define("leadJet_medConstituentHcalDepthEF2", "goodJet_medConstituentHcalDepthEF2[0]")
    df = df.Define("leadJet_medConstituentHcalDepthEF3", "goodJet_medConstituentHcalDepthEF3[0]")
    df = df.Define("leadJet_medConstituentHcalDepthEF4", "goodJet_medConstituentHcalDepthEF4[0]")
    df = df.Define("leadJet_medConstituentHcalDepthEF5", "goodJet_medConstituentHcalDepthEF5[0]")
    df = df.Define("leadJet_medConstituentHcalDepthEF6", "goodJet_medConstituentHcalDepthEF6[0]")
    df = df.Define("leadJet_medConstituentHcalDepthEF7", "goodJet_medConstituentHcalDepthEF7[0]")
    
    df = df.Define("leadJet_minConstituentHcalDepthEF1", "goodJet_minConstituentHcalDepthEF1[0]")
    df = df.Define("leadJet_minConstituentHcalDepthEF2", "goodJet_minConstituentHcalDepthEF2[0]")
    df = df.Define("leadJet_minConstituentHcalDepthEF3", "goodJet_minConstituentHcalDepthEF3[0]")
    df = df.Define("leadJet_minConstituentHcalDepthEF4", "goodJet_minConstituentHcalDepthEF4[0]")
    df = df.Define("leadJet_minConstituentHcalDepthEF5", "goodJet_minConstituentHcalDepthEF5[0]")
    df = df.Define("leadJet_minConstituentHcalDepthEF6", "goodJet_minConstituentHcalDepthEF6[0]")
    df = df.Define("leadJet_minConstituentHcalDepthEF7", "goodJet_minConstituentHcalDepthEF7[0]")
    
    df = df.Define("leadJet_maxConstituentHcalDepthEF1", "goodJet_maxConstituentHcalDepthEF1[0]")
    df = df.Define("leadJet_maxConstituentHcalDepthEF2", "goodJet_maxConstituentHcalDepthEF2[0]")
    df = df.Define("leadJet_maxConstituentHcalDepthEF3", "goodJet_maxConstituentHcalDepthEF3[0]")
    df = df.Define("leadJet_maxConstituentHcalDepthEF4", "goodJet_maxConstituentHcalDepthEF4[0]")
    df = df.Define("leadJet_maxConstituentHcalDepthEF5", "goodJet_maxConstituentHcalDepthEF5[0]")
    df = df.Define("leadJet_maxConstituentHcalDepthEF6", "goodJet_maxConstituentHcalDepthEF6[0]")
    df = df.Define("leadJet_maxConstituentHcalDepthEF7", "goodJet_maxConstituentHcalDepthEF7[0]")
    
    # separate values into barrel and endcap regions
    df = df.Define("leadJet_HB_mask", "abs(leadJet_Eta)<1.305")
    df = df.Define("leadJet_HE_mask", "(abs(leadJet_Eta)<2.5) && (abs(leadJet_Eta)>1.305)")
    dfHB = df.Filter("leadJet_HB_mask")
    dfHE = df.Filter("leadJet_HE_mask")

    return df, dfHB, dfHE

def createHisto1D(
    hist_list,
    df,
    hist_name,
    var, n_bins, bin_low, bin_high,
):

    histoModel = ROOT.RDF.TH1DModel(
        name=hist_name,
        title=hist_name,
        nbinsx=n_bins,
        xlow=bin_low,
        xup=bin_high
    )

    hist_list.append(df.Histo1D(histoModel, var))

    return hist_list

def createHisto2D(
    hist_list,
    df,
    hist_name,
    var_x, n_bins_x, bin_low_x, bin_high_x,
    var_y, n_bins_y, bin_low_y, bin_high_y,
):

    histoModel = ROOT.RDF.TH2DModel(
        name=hist_name,
        title=hist_name,
        nbinsx=n_bins_x,
        xlow=bin_low_x,
        xup=bin_high_x,
        nbinsy=n_bins_y,
        ylow=bin_low_y,
        yup=bin_high_y,
    )

    hist_list.append(df.Histo2D(histoModel, var_x, var_y))

    return hist_list


def bookHistos(
    df, dfHB, dfHE
):

    hists = []

    hists = createHisto1D(hists, df, "leadJet_Pt", "leadJet_Pt", 200, 0, 1000)
    hists = createHisto1D(hists, df, "leadJet_Eta", "leadJet_Eta", 100, -3, 3)
    hists = createHisto1D(hists, df, "leadJet_Phi", "leadJet_Phi", 100, -3.14, 3.14)
    hists = createHisto1D(hists, df, "leadJet_E", "leadJet_E", 200, 0, 1000)
    hists = createHisto1D(hists, df, "leadJet_nConstituent", "leadJet_nConstituent", 200, 0, 200)
    
    hists = createHisto2D(hists, dfHB, "leadJet_HB_avgConstituentHcalDepthEF1", "leadJet_avgConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHB, "leadJet_HB_avgConstituentHcalDepthEF2", "leadJet_avgConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHB, "leadJet_HB_avgConstituentHcalDepthEF3", "leadJet_avgConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHB, "leadJet_HB_avgConstituentHcalDepthEF4", "leadJet_avgConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_avgConstituentHcalDepthEF1", "leadJet_avgConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_avgConstituentHcalDepthEF2", "leadJet_avgConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_avgConstituentHcalDepthEF3", "leadJet_avgConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_avgConstituentHcalDepthEF4", "leadJet_avgConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_avgConstituentHcalDepthEF5", "leadJet_avgConstituentHcalDepthEF5", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_avgConstituentHcalDepthEF6", "leadJet_avgConstituentHcalDepthEF6", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_avgConstituentHcalDepthEF7", "leadJet_avgConstituentHcalDepthEF7", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    
    hists = createHisto2D(hists, dfHB, "leadJet_HB_pTWeightedAvgConstituentHcalDepthEF1", "leadJet_pTWeightedAvgConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHB, "leadJet_HB_pTWeightedAvgConstituentHcalDepthEF2", "leadJet_pTWeightedAvgConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHB, "leadJet_HB_pTWeightedAvgConstituentHcalDepthEF3", "leadJet_pTWeightedAvgConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHB, "leadJet_HB_pTWeightedAvgConstituentHcalDepthEF4", "leadJet_pTWeightedAvgConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_pTWeightedAvgConstituentHcalDepthEF1", "leadJet_pTWeightedAvgConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_pTWeightedAvgConstituentHcalDepthEF2", "leadJet_pTWeightedAvgConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_pTWeightedAvgConstituentHcalDepthEF3", "leadJet_pTWeightedAvgConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_pTWeightedAvgConstituentHcalDepthEF4", "leadJet_pTWeightedAvgConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_pTWeightedAvgConstituentHcalDepthEF5", "leadJet_pTWeightedAvgConstituentHcalDepthEF5", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_pTWeightedAvgConstituentHcalDepthEF6", "leadJet_pTWeightedAvgConstituentHcalDepthEF6", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_pTWeightedAvgConstituentHcalDepthEF7", "leadJet_pTWeightedAvgConstituentHcalDepthEF7", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    
    hists = createHisto2D(hists, dfHB, "leadJet_HB_medConstituentHcalDepthEF1", "leadJet_medConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHB, "leadJet_HB_medConstituentHcalDepthEF2", "leadJet_medConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHB, "leadJet_HB_medConstituentHcalDepthEF3", "leadJet_medConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHB, "leadJet_HB_medConstituentHcalDepthEF4", "leadJet_medConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_medConstituentHcalDepthEF1", "leadJet_medConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_medConstituentHcalDepthEF2", "leadJet_medConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_medConstituentHcalDepthEF3", "leadJet_medConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_medConstituentHcalDepthEF4", "leadJet_medConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_medConstituentHcalDepthEF5", "leadJet_medConstituentHcalDepthEF5", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_medConstituentHcalDepthEF6", "leadJet_medConstituentHcalDepthEF6", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_medConstituentHcalDepthEF7", "leadJet_medConstituentHcalDepthEF7", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    
    hists = createHisto2D(hists, dfHB, "leadJet_HB_minConstituentHcalDepthEF1", "leadJet_minConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHB, "leadJet_HB_minConstituentHcalDepthEF2", "leadJet_minConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHB, "leadJet_HB_minConstituentHcalDepthEF3", "leadJet_minConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHB, "leadJet_HB_minConstituentHcalDepthEF4", "leadJet_minConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_minConstituentHcalDepthEF1", "leadJet_minConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_minConstituentHcalDepthEF2", "leadJet_minConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_minConstituentHcalDepthEF3", "leadJet_minConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_minConstituentHcalDepthEF4", "leadJet_minConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_minConstituentHcalDepthEF5", "leadJet_minConstituentHcalDepthEF5", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_minConstituentHcalDepthEF6", "leadJet_minConstituentHcalDepthEF6", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_minConstituentHcalDepthEF7", "leadJet_minConstituentHcalDepthEF7", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    
    hists = createHisto2D(hists, dfHB, "leadJet_HB_maxConstituentHcalDepthEF1", "leadJet_maxConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHB, "leadJet_HB_maxConstituentHcalDepthEF2", "leadJet_maxConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHB, "leadJet_HB_maxConstituentHcalDepthEF3", "leadJet_maxConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHB, "leadJet_HB_maxConstituentHcalDepthEF4", "leadJet_maxConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_maxConstituentHcalDepthEF1", "leadJet_maxConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_maxConstituentHcalDepthEF2", "leadJet_maxConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_maxConstituentHcalDepthEF3", "leadJet_maxConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_maxConstituentHcalDepthEF4", "leadJet_maxConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_maxConstituentHcalDepthEF5", "leadJet_maxConstituentHcalDepthEF5", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_maxConstituentHcalDepthEF6", "leadJet_maxConstituentHcalDepthEF6", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, dfHE, "leadJet_HE_maxConstituentHcalDepthEF7", "leadJet_maxConstituentHcalDepthEF7", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)

    return hists



def runDistributed(
    args,
    inputChain,
):

    with Client(args.scheduler_address) as client:
        df = ROOT.RDF.Experimental.Distributed.Dask.RDataFrame(
            inputChain,
            daskclient=client,
            npartitions=args.npartitions,
        )

        histoList = bookHistos(df)

        ROOT.RDF.Experimental.Distributed.RunGraphs(histoList)

    return histoList

def runLocal(
    args,
    inputChain,
):

    print("Running locally")
    
    declareCppMethods()
    print("Declared Cpp Methods using ROOT gInterpreter")
    
    df = ROOT.RDataFrame(inputChain)
    print("Created dataframe")
    print("number of events: ", df.Count().GetValue())
    
    df, dfHB, dfHE = prepareDataframe(df)
    print("Prepared columns in dataframe")

    histoList = bookHistos(df, dfHB, dfHE)
    print("Booked histos")

    return histoList

def postProcessResults(results):

    print("Post-processing")
    newResults = []
    for i, res in enumerate(results):
        print(f"    processing hist #{i}")
        if hasattr(res, "GetValue"):
            h = res.GetValue()
            newResults.append(h)
        else:
            assert hasattr(res, "GetKeys")
            for key in res.GetKeys():
                h = res[key]
                newResults.append(h)
    return newResults

def saveHistos(results, outputFileName):
    print("Saving to ROOT file")
    with ROOT.TFile.Open(outputFileName, "recreate") as outputFile:
        for result in results:
            outputFile.WriteObject(result, result.GetName())



###################################################################################################
# MAIN METHOD

def main():

    startTime = time()
    args = parseArgs()

    ROOT.TH1.AddDirectory(False)
    ROOT.gROOT.SetBatch(True)

    print("Loading json")
    with open(args.input) as f:
        fileDict = json.load(f)

    filesForDF = fileDict[args.dataset]

    print("Creating TChain")
    inputChain = createChain(filesForDF, treesForDF)

    # results = runDistributed(args, inputChain)
    results = runLocal(args, inputChain)

    results = postProcessResults(results)

    saveHistos(results, args.output)

    print(f"Histograms created and saved in {time() - startTime:.2f} seconds")


###################################################################################################
# RUN SCRIPT

if __name__ == "__main__":
    main()
