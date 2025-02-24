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

    # function for unflattening jet constituent info
    ROOT.gInterpreter.Declare("""
#include <vector>

std::vector<std::vector<double>> reconstructJetConstituentInfo(const ROOT::VecOps::RVec<double>& flat_data, const ROOT::VecOps::RVec<double>& sizes) {
    std::vector<std::vector<double>> perJetData;
    int index = 0;
    
    for (int i = 0; i < sizes.size(); ++i) {
        int size = (int)sizes[i];
        std::vector<double> jetConstituents(flat_data.begin() + index, flat_data.begin() + index + size);
        perJetData.push_back(jetConstituents);
        index += size;
    }
    
    return perJetData;
}
""")

    # function for applying jet mask to nested variable
    ROOT.gInterpreter.Declare(R"""
#include <vector>

std::vector<std::vector<double>> filterJetsByMask(const std::vector<std::vector<double>>& jetData, const ROOT::VecOps::RVec<int>& mask) {
    std::vector<std::vector<double>> filteredJets;
    for (size_t i = 0; i < mask.size(); ++i) {
        if (mask[i]) {
            filteredJets.push_back(jetData[i]);
        }
    }
    return filteredJets;
}
""")

    # function for computing averages of constituent info per jet
    ROOT.gInterpreter.Declare("""
#include <vector>

std::vector<double> computeJetAverages(const std::vector<std::vector<double>>& perJetData) {
    std::vector<double> jetAverages;
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

std::vector<std::vector<double>> computeJetConstituentPt(
    const std::vector<std::vector<double>>& JetConstituentPx,
    const std::vector<std::vector<double>>& JetConstituentPy)
{
    std::vector<std::vector<double>> JetConstituentPt;
    for (size_t i = 0; i < JetConstituentPx.size(); ++i) {
        std::vector<double> ptValues;
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

std::vector<double> computePtWeightedJetAverages(
    const std::vector<std::vector<double>>& perJetConstituentVar,
    const std::vector<std::vector<double>>& perJetConstituentPt)
{
    std::vector<double> jetWeightedAverages;
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
    
std::vector<double> computeJetMedians(const std::vector<std::vector<double>>& perJetConstituentVar) {
    std::vector<double> jetMedians;
    for(size_t i = 0; i < perJetConstituentVar.size(); ++i) {
        if (!perJetConstituentVar[i].empty()) {
            std::vector<double> sortedData = perJetConstituentVar[i];
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
    
std::vector<double> computeJetMinimums(const std::vector<std::vector<double>>& perJetConstituentVar) {
    std::vector<double> jetMinimums;
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
    
std::vector<double> computeJetMaximums(const std::vector<std::vector<double>>& perJetConstituentVar) {
    std::vector<double> jetMaximums;
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

    return

def prepareDataframe(df):

    # include events with at least 1 jet with pT > 20 GeV
    df = df.Define("Jet_mask", "(Jet_Pt > 20) && (abs(Jet_Eta)<2.5)") # apply jet selections here
    df = df.Filter("Sum(Jet_mask)>=1") # require at least one good jet per event
    print("Filtered events to contain at least one good jet")

    # unflatten jet constituent info
    df = df.Define("Jet_Constituent_Px", "reconstructJetConstituentInfo(Constituent_Px, Jet_nConstituent)")
    df = df.Define("Jet_Constituent_Py", "reconstructJetConstituentInfo(Constituent_Py, Jet_nConstituent)")
    df = df.Define("Jet_Constituent_Pz", "reconstructJetConstituentInfo(Constituent_Pz, Jet_nConstituent)")
    df = df.Define("Jet_Constituent_E", "reconstructJetConstituentInfo(Constituent_E, Jet_nConstituent)")
    print("Unflatted jet constituent kinematics variables")
    
    df = df.Define("Jet_Constituent_hcalDepthEF1", "reconstructJetConstituentInfo(Constituent_hcalDepthEF1, Jet_nConstituent)")
    df = df.Define("Jet_Constituent_hcalDepthEF2", "reconstructJetConstituentInfo(Constituent_hcalDepthEF2, Jet_nConstituent)")
    df = df.Define("Jet_Constituent_hcalDepthEF3", "reconstructJetConstituentInfo(Constituent_hcalDepthEF3, Jet_nConstituent)")
    df = df.Define("Jet_Constituent_hcalDepthEF4", "reconstructJetConstituentInfo(Constituent_hcalDepthEF4, Jet_nConstituent)")
    df = df.Define("Jet_Constituent_hcalDepthEF5", "reconstructJetConstituentInfo(Constituent_hcalDepthEF5, Jet_nConstituent)")
    df = df.Define("Jet_Constituent_hcalDepthEF6", "reconstructJetConstituentInfo(Constituent_hcalDepthEF6, Jet_nConstituent)")
    df = df.Define("Jet_Constituent_hcalDepthEF7", "reconstructJetConstituentInfo(Constituent_hcalDepthEF7, Jet_nConstituent)")
    print("Unflatted jet constituent HCAL depth fractions")
    
    # filter out bad jets
    df = df.Define("goodJet_Pt", "Jet_Pt[Jet_mask]")
    df = df.Define("goodJet_Eta", "Jet_Eta[Jet_mask]")
    df = df.Define("goodJet_Phi", "Jet_Phi[Jet_mask]")
    df = df.Define("goodJet_E", "Jet_E[Jet_mask]")
    df = df.Define("goodJet_nConstituent", "Jet_nConstituent[Jet_mask]")
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
    df = df.Define("leadJet_nConstituent", "goodJet_nConstituent[0]")
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
    
    

    return df

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
    df
):

    hists = []

    hists = createHisto1D(hists, df, "leadJet_Pt", "leadJet_Pt", 200, 0, 1000)
    hists = createHisto1D(hists, df, "leadJet_Eta", "leadJet_Eta", 100, -3, 3)
    hists = createHisto1D(hists, df, "leadJet_Phi", "leadJet_Phi", 100, -3.14, 3.14)
    hists = createHisto1D(hists, df, "leadJet_E", "leadJet_E", 200, 0, 1000)
    hists = createHisto1D(hists, df, "leadJet_nConstituent", "leadJet_nConstituent", 200, 0, 200)
    hists = createHisto2D(hists, df, "leadJet_avgConstituentHcalDepthEF1", "leadJet_avgConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_avgConstituentHcalDepthEF2", "leadJet_avgConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_avgConstituentHcalDepthEF3", "leadJet_avgConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_avgConstituentHcalDepthEF4", "leadJet_avgConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_avgConstituentHcalDepthEF5", "leadJet_avgConstituentHcalDepthEF5", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_avgConstituentHcalDepthEF6", "leadJet_avgConstituentHcalDepthEF6", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_avgConstituentHcalDepthEF7", "leadJet_avgConstituentHcalDepthEF7", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_pTWeightedAvgConstituentHcalDepthEF1", "leadJet_pTWeightedAvgConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_pTWeightedAvgConstituentHcalDepthEF2", "leadJet_pTWeightedAvgConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_pTWeightedAvgConstituentHcalDepthEF3", "leadJet_pTWeightedAvgConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_pTWeightedAvgConstituentHcalDepthEF4", "leadJet_pTWeightedAvgConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_pTWeightedAvgConstituentHcalDepthEF5", "leadJet_pTWeightedAvgConstituentHcalDepthEF5", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_pTWeightedAvgConstituentHcalDepthEF6", "leadJet_pTWeightedAvgConstituentHcalDepthEF6", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_pTWeightedAvgConstituentHcalDepthEF7", "leadJet_pTWeightedAvgConstituentHcalDepthEF7", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_medConstituentHcalDepthEF1", "leadJet_medConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_medConstituentHcalDepthEF2", "leadJet_medConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_medConstituentHcalDepthEF3", "leadJet_medConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_medConstituentHcalDepthEF4", "leadJet_medConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_medConstituentHcalDepthEF5", "leadJet_medConstituentHcalDepthEF5", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_medConstituentHcalDepthEF6", "leadJet_medConstituentHcalDepthEF6", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_medConstituentHcalDepthEF7", "leadJet_medConstituentHcalDepthEF7", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_medConstituentHcalDepthEF1", "leadJet_medConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_medConstituentHcalDepthEF2", "leadJet_medConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_medConstituentHcalDepthEF3", "leadJet_medConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_medConstituentHcalDepthEF4", "leadJet_medConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_medConstituentHcalDepthEF5", "leadJet_medConstituentHcalDepthEF5", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_medConstituentHcalDepthEF6", "leadJet_medConstituentHcalDepthEF6", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_medConstituentHcalDepthEF7", "leadJet_medConstituentHcalDepthEF7", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_minConstituentHcalDepthEF1", "leadJet_minConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_minConstituentHcalDepthEF2", "leadJet_minConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_minConstituentHcalDepthEF3", "leadJet_minConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_minConstituentHcalDepthEF4", "leadJet_minConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_minConstituentHcalDepthEF5", "leadJet_minConstituentHcalDepthEF5", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_minConstituentHcalDepthEF6", "leadJet_minConstituentHcalDepthEF6", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_minConstituentHcalDepthEF7", "leadJet_minConstituentHcalDepthEF7", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_maxConstituentHcalDepthEF1", "leadJet_maxConstituentHcalDepthEF1", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_maxConstituentHcalDepthEF2", "leadJet_maxConstituentHcalDepthEF2", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_maxConstituentHcalDepthEF3", "leadJet_maxConstituentHcalDepthEF3", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_maxConstituentHcalDepthEF4", "leadJet_maxConstituentHcalDepthEF4", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_maxConstituentHcalDepthEF5", "leadJet_maxConstituentHcalDepthEF5", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_maxConstituentHcalDepthEF6", "leadJet_maxConstituentHcalDepthEF6", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
    hists = createHisto2D(hists, df, "leadJet_maxConstituentHcalDepthEF7", "leadJet_maxConstituentHcalDepthEF7", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)

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
    
    df = prepareDataframe(df)
    print("Prepared columns in dataframe")

    histoList = bookHistos(df)
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
