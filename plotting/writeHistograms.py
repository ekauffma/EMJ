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
    "EMJMCNtuplizer/Events",
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

RVec<double> computeJet_avg(const RVec<RVec<double>>& perJetData) {
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

RVec<double> computeJet_pTWeightedAvg(
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
    
RVec<double> computeJet_med(const RVec<RVec<double>>& perJetConstituentVar) {
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
    
RVec<double> computeJet_min(const RVec<RVec<double>>& perJetConstituentVar) {
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
    
RVec<double> computeJet_max(const RVec<RVec<double>>& perJetConstituentVar) {
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

    # function to compute event-level median
    ROOT.gInterpreter.Declare("""
template <typename T> 
T Median(ROOT::VecOps::RVec<T> v) {
    if (v.empty()) return T();
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    if (n % 2 == 0)
        return (v[n/2 - 1] + v[n/2]) / T(2);
    else
        return v[n/2];
}
""")
    
    # function to compute event-level fraction
    ROOT.gInterpreter.Declare("""
template <typename T>
double EventLevelFractionEqual(ROOT::VecOps::RVec<T> v, T val) {
    if (v.empty()) return -1.0;
    int count = 0;
    for (const auto& x : v) {
        if (x==val) ++count;
    }
    return static_cast<double>(count) / v.size();
}
""")

    ROOT.gInterpreter.Declare("""
template <typename T, typename U>
double EventLevelEnergyWeightedFractionEqual(
    const ROOT::VecOps::RVec<T>& v,
    const ROOT::VecOps::RVec<U>& energies,
    T val
) {
    if (v.size() != energies.size() || v.empty())
        return -1.0;

    double total_energy = 0.0;
    double matched_energy = 0.0;

    for (size_t i = 0; i < v.size(); ++i) {
        total_energy += energies[i];
        if (v[i] == val) {
            matched_energy += energies[i];
        }
    }

    if (total_energy > 0.0)
        return matched_energy / total_energy;
    else
        return 0.0;
}
""")


    ROOT.gInterpreter.Declare("""
template <typename T>
ROOT::VecOps::RVec<double> JetLevelFractionEqual(ROOT::VecOps::RVec<ROOT::VecOps::RVec<T>> jets, T val) {
    ROOT::VecOps::RVec<double> fractions;
    for (const auto& jet : jets) {
        if (jet.empty()) {
            fractions.push_back(-1.0);
            continue;
        }
        int count = 0;
        for (const auto& x : jet) {
            if (x == val) ++count;
        }
        fractions.push_back(static_cast<double>(count) / jet.size());
    }
    return fractions;
}
""")

    ROOT.gInterpreter.Declare("""
template <typename T, typename U>
ROOT::VecOps::RVec<double> JetLevelEnergyWeightedFractionEqual(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<T>>& jets,
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<U>>& energies,
    const T& val
) {
    ROOT::VecOps::RVec<double> result;
    size_t n_jets = jets.size();

    for (size_t j = 0; j < n_jets; ++j) {
        const auto& jet = jets[j];
        const auto& energy = energies[j];

        if (jet.size() != energy.size() || jet.empty()) {
            result.push_back(-1.0);
            continue;
        }

        double total_energy = 0.0;
        double matched_energy = 0.0;

        for (size_t i = 0; i < jet.size(); ++i) {
            total_energy += energy[i];
            if (jet[i] == val) {
                matched_energy += energy[i];
            }
        }

        result.push_back(total_energy > 0.0 ? matched_energy / total_energy : 0.0);
    }

    return result;
}
""")

    ROOT.gInterpreter.Declare("""
#include <vector>
#include <cmath>
using namespace ROOT::VecOps;
    
template <typename T>
RVec<RVec<T>> reconstructNestedArray(const RVec<T>& vec, const RVec<int>& counts) {
    RVec<RVec<T>> nested;
    int offset = 0;
    for (size_t i = 0; i < counts.size(); ++i) {
        int n = counts[i];
        nested.emplace_back(vec.begin() + offset, vec.begin() + offset + n);
        offset += n;
    }
    return nested;
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
    for i in range(1,8):
        df = df.Define(f"filteredConstituent_hcalDepthEF{i}", f"filterVector(Constituent_hcalDepthEF{i}, valid_indices)")
    df = df.Define("filteredConstituent_rawHoverE", "filterVector(Constituent_rawHoverE, valid_indices)")
    df = df.Define("filteredConstituent_HoverE", "filterVector(Constituent_HoverE, valid_indices)")
    print("Filtered out constituents with no deposits in HCAL")
    
    
    # unflatten jet constituent info
    df = df.Define("Jet_Constituent_Px", "reconstructNestedArray(filteredConstituent_Px, Jet_nFilteredConstituent)")
    df = df.Define("Jet_Constituent_Py", "reconstructNestedArray(filteredConstituent_Py, Jet_nFilteredConstituent)")
    df = df.Define("Jet_Constituent_Pz", "reconstructNestedArray(filteredConstituent_Pz, Jet_nFilteredConstituent)")
    df = df.Define("Jet_Constituent_E", "reconstructNestedArray(filteredConstituent_E, Jet_nFilteredConstituent)")
    print("Unflatted jet constituent kinematics variables")
    
    for i in range(1,8):
        df = df.Define(f"Jet_Constituent_hcalDepthEF{i}", f"reconstructNestedArray(filteredConstituent_hcalDepthEF{i}, Jet_nFilteredConstituent)")
    df = df.Define("Jet_Constituent_rawHoverE", "reconstructNestedArray(filteredConstituent_rawHoverE, Jet_nFilteredConstituent)")
    df = df.Define("Jet_Constituent_HoverE", "reconstructNestedArray(filteredConstituent_HoverE, Jet_nFilteredConstituent)")
    print("Unflatted jet constituent HCAL depth fractions")
    
    # unflatten jet rechit info
    df = df.Define("Jet_JetRecHit_auxTDC3", "reconstructNestedArray(JetRecHit_auxTDC3, Jet_nJetRecHit)")
    df = df.Define("Jet_JetRecHit_energy", "reconstructNestedArray(JetRecHit_energy, Jet_nJetRecHit)")
    df = df.Define("Jet_JetRecHit_eta", "reconstructNestedArray(JetRecHit_eta, Jet_nJetRecHit)")
    df = df.Define("Jet_JetRecHit_phi", "reconstructNestedArray(JetRecHit_phi, Jet_nJetRecHit)")
    
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
    for i in range(1,8):
        df = df.Define(f"goodJet_Constituent_hcalDepthEF{i}", f"filterJetsByMask(Jet_Constituent_hcalDepthEF{i}, Jet_mask)")
    df = df.Define("goodJet_Constituent_rawHoverE", "filterJetsByMask(Jet_Constituent_rawHoverE, Jet_mask)")
    df = df.Define("goodJet_Constituent_HoverE", "filterJetsByMask(Jet_Constituent_HoverE, Jet_mask)")
    df = df.Define("goodJet_JetRecHit_auxTDC3", "filterJetsByMask(Jet_JetRecHit_auxTDC3, Jet_mask)")
    df = df.Define("goodJet_JetRecHit_energy", "filterJetsByMask(Jet_JetRecHit_energy, Jet_mask)")
    df = df.Define("goodJet_JetRecHit_eta", "filterJetsByMask(Jet_JetRecHit_eta, Jet_mask)")
    df = df.Define("goodJet_JetRecHit_phi", "filterJetsByMask(Jet_JetRecHit_phi, Jet_mask)")
    print("Filtered out bad jets")
    
    # compute jet constituent Pt
    df = df.Define("goodJet_Constituent_Pt", "computeJetConstituentPt(goodJet_Constituent_Px, goodJet_Constituent_Py)")
    print("Defined jet constituent transverse momentum")
    
    for metric in ["avg", "pTWeightedAvg", "med", "min", "max"]:
        for i in range(1,8):
            if metric=="pTWeightedAvg":
                df = df.Define(f"goodJet_{metric}ConstituentHcalDepthEF{i}", f"computeJet_{metric}(goodJet_Constituent_hcalDepthEF{i}, goodJet_Constituent_Pt)")
            else:
                df = df.Define(f"goodJet_{metric}ConstituentHcalDepthEF{i}", f"computeJet_{metric}(goodJet_Constituent_hcalDepthEF{i})")
        if metric=="pTWeightedAvg":
            df = df.Define(f"goodJet_{metric}ConstituentRawHoverE", f"computeJet_{metric}(goodJet_Constituent_rawHoverE, goodJet_Constituent_Pt)")
            df = df.Define(f"goodJet_{metric}ConstituentHoverE", f"computeJet_{metric}(goodJet_Constituent_HoverE, goodJet_Constituent_Pt)")
        else:
            df = df.Define(f"goodJet_{metric}ConstituentRawHoverE", f"computeJet_{metric}(goodJet_Constituent_rawHoverE)")
            df = df.Define(f"goodJet_{metric}ConstituentHoverE", f"computeJet_{metric}(goodJet_Constituent_HoverE)")
    df = df.Redefine("goodJet_JetRecHit_auxTDC3", "ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>(goodJet_JetRecHit_auxTDC3)")
    df = df.Define("goodJet_JetRecHit_auxTDC3_frac0", "JetLevelFractionEqual(goodJet_JetRecHit_auxTDC3, 0)")
    df = df.Define("goodJet_JetRecHit_auxTDC3_frac1", "JetLevelFractionEqual(goodJet_JetRecHit_auxTDC3, 1)")
    df = df.Define("goodJet_JetRecHit_auxTDC3_frac2", "JetLevelFractionEqual(goodJet_JetRecHit_auxTDC3, 2)")
    df = df.Define("goodJet_JetRecHit_auxTDC3_frac3", "JetLevelFractionEqual(goodJet_JetRecHit_auxTDC3, 3)")
    df = df.Define("goodJet_JetRecHit_auxTDC3_energyWeightedFrac0", "JetLevelEnergyWeightedFractionEqual(goodJet_JetRecHit_auxTDC3, goodJet_JetRecHit_energy, 0)")
    df = df.Define("goodJet_JetRecHit_auxTDC3_energyWeightedFrac1", "JetLevelEnergyWeightedFractionEqual(goodJet_JetRecHit_auxTDC3, goodJet_JetRecHit_energy, 1)")
    df = df.Define("goodJet_JetRecHit_auxTDC3_energyWeightedFrac2", "JetLevelEnergyWeightedFractionEqual(goodJet_JetRecHit_auxTDC3, goodJet_JetRecHit_energy, 2)")
    df = df.Define("goodJet_JetRecHit_auxTDC3_energyWeightedFrac3", "JetLevelEnergyWeightedFractionEqual(goodJet_JetRecHit_auxTDC3, goodJet_JetRecHit_energy, 3)")
    df = df.Define("goodJet_JetRecHit_energy_avg", "computeJet_avg(goodJet_JetRecHit_energy)")
    df = df.Define("goodJet_JetRecHit_energy_max", "computeJet_max(goodJet_JetRecHit_energy)")
    df = df.Define("goodJet_JetRecHit_energy_min", "computeJet_min(goodJet_JetRecHit_energy)")
    
    # define leading jet values
    df = df.Define("leadJet_Pt", "goodJet_Pt[0]")
    df = df.Define("leadJet_Eta", "goodJet_Eta[0]")
    df = df.Define("leadJet_Phi", "goodJet_Phi[0]")
    df = df.Define("leadJet_E", "goodJet_E[0]")
    df = df.Define("leadJet_nConstituent", "goodJet_nFilteredConstituent[0]")
    
    # filter out events where the leading jet has no remaining constituents
    df = df.Filter("leadJet_nConstituent > 0")
    
    for metric in ["avg", "pTWeightedAvg", "med", "min", "max"]:
        for i in range(1,8):
            df = df.Define(f"leadJet_{metric}ConstituentHcalDepthEF{i}", f"goodJet_{metric}ConstituentHcalDepthEF{i}[0]")
        df = df.Define(f"leadJet_{metric}ConstituentRawHoverE", f"goodJet_{metric}ConstituentRawHoverE[0]")
        df = df.Define(f"leadJet_{metric}ConstituentHoverE", f"goodJet_{metric}ConstituentHoverE[0]")
    df = df.Define(f"leadJet_JetRecHit_auxTDC3_frac0", f"goodJet_JetRecHit_auxTDC3_frac0[0]")
    df = df.Define(f"leadJet_JetRecHit_auxTDC3_frac1", f"goodJet_JetRecHit_auxTDC3_frac1[0]")
    df = df.Define(f"leadJet_JetRecHit_auxTDC3_frac2", f"goodJet_JetRecHit_auxTDC3_frac2[0]")
    df = df.Define(f"leadJet_JetRecHit_auxTDC3_frac3", f"goodJet_JetRecHit_auxTDC3_frac3[0]")
    df = df.Define(f"leadJet_JetRecHit_auxTDC3_energyWeightedFrac0", f"goodJet_JetRecHit_auxTDC3_energyWeightedFrac0[0]")
    df = df.Define(f"leadJet_JetRecHit_auxTDC3_energyWeightedFrac1", f"goodJet_JetRecHit_auxTDC3_energyWeightedFrac1[0]")
    df = df.Define(f"leadJet_JetRecHit_auxTDC3_energyWeightedFrac2", f"goodJet_JetRecHit_auxTDC3_energyWeightedFrac2[0]")
    df = df.Define(f"leadJet_JetRecHit_auxTDC3_energyWeightedFrac3", f"goodJet_JetRecHit_auxTDC3_energyWeightedFrac3[0]")
    df = df.Define(f"leadJet_JetRecHit_energy_avg", f"goodJet_JetRecHit_energy_avg[0]")
    df = df.Define(f"leadJet_JetRecHit_energy_max", f"goodJet_JetRecHit_energy_max[0]")
    df = df.Define(f"leadJet_JetRecHit_energy_min", f"goodJet_JetRecHit_energy_min[0]")

    df = df.Define("RecHit_E_max", "ROOT::VecOps::Max(RecHit_E)")
    df = df.Define("RecHit_E_min", "ROOT::VecOps::Min(RecHit_E)")
    df = df.Define("RecHit_E_median", "Median(RecHit_E)")
    df = df.Define("RecHit_E_mean", "ROOT::VecOps::Mean(RecHit_E)")
    df = df.Define("RecHit_auxTDC3_frac0", "EventLevelFractionEqual(RecHit_auxTDC3, 0)")
    df = df.Define("RecHit_auxTDC3_frac1", "EventLevelFractionEqual(RecHit_auxTDC3, 1)")
    df = df.Define("RecHit_auxTDC3_frac2", "EventLevelFractionEqual(RecHit_auxTDC3, 2)")
    df = df.Define("RecHit_auxTDC3_frac3", "EventLevelFractionEqual(RecHit_auxTDC3, 3)")
    df = df.Define("RecHit_auxTDC3_energyWeightedFrac0", "EventLevelEnergyWeightedFractionEqual(RecHit_auxTDC3, RecHit_E, 0)")
    df = df.Define("RecHit_auxTDC3_energyWeightedFrac1", "EventLevelEnergyWeightedFractionEqual(RecHit_auxTDC3, RecHit_E, 1)")
    df = df.Define("RecHit_auxTDC3_energyWeightedFrac2", "EventLevelEnergyWeightedFractionEqual(RecHit_auxTDC3, RecHit_E, 2)")
    df = df.Define("RecHit_auxTDC3_energyWeightedFrac3", "EventLevelEnergyWeightedFractionEqual(RecHit_auxTDC3, RecHit_E, 3)")

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

    hists = createHisto1D(hists, df, "nPV", "nPV", 100, 0, 100)
    hists = createHisto1D(hists, df, "leadJet_Pt", "leadJet_Pt", 200, 0, 1000)
    hists = createHisto1D(hists, df, "leadJet_Eta", "leadJet_Eta", 100, -3, 3)
    hists = createHisto1D(hists, df, "leadJet_Phi", "leadJet_Phi", 100, -3.14, 3.14)
    hists = createHisto1D(hists, df, "leadJet_E", "leadJet_E", 200, 0, 1000)
    hists = createHisto1D(hists, df, "leadJet_nConstituent", "leadJet_nConstituent", 200, 0, 200)
    hists = createHisto1D(hists, df, "leadJet_JetRecHit_auxTDC3_frac0", "leadJet_JetRecHit_auxTDC3_frac0", 200, -1, 1)
    hists = createHisto1D(hists, df, "leadJet_JetRecHit_auxTDC3_frac1", "leadJet_JetRecHit_auxTDC3_frac1", 200, -1, 1)
    hists = createHisto1D(hists, df, "leadJet_JetRecHit_auxTDC3_frac2", "leadJet_JetRecHit_auxTDC3_frac2", 200, -1, 1)
    hists = createHisto1D(hists, df, "leadJet_JetRecHit_auxTDC3_frac3", "leadJet_JetRecHit_auxTDC3_frac3", 200, -1, 1)
    hists = createHisto1D(hists, df, "leadJet_JetRecHit_auxTDC3_energyWeightedFrac0", "leadJet_JetRecHit_auxTDC3_energyWeightedFrac0", 200, -1, 1)
    hists = createHisto1D(hists, df, "leadJet_JetRecHit_auxTDC3_energyWeightedFrac1", "leadJet_JetRecHit_auxTDC3_energyWeightedFrac1", 200, -1, 1)
    hists = createHisto1D(hists, df, "leadJet_JetRecHit_auxTDC3_energyWeightedFrac2", "leadJet_JetRecHit_auxTDC3_energyWeightedFrac2", 200, -1, 1)
    hists = createHisto1D(hists, df, "leadJet_JetRecHit_auxTDC3_energyWeightedFrac3", "leadJet_JetRecHit_auxTDC3_energyWeightedFrac3", 200, -1, 1)
    hists = createHisto1D(hists, df, "leadJet_JetRecHit_energy_avg", "leadJet_JetRecHit_energy_avg", 100, 0, 4)
    hists = createHisto1D(hists, df, "leadJet_JetRecHit_energy_max", "leadJet_JetRecHit_energy_max", 400, 0, 200)
    hists = createHisto1D(hists, df, "leadJet_JetRecHit_energy_min", "leadJet_JetRecHit_energy_min", 100, 0, 2)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_Pt", "leadJet_Pt", 200, 0, 1000)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_Eta", "leadJet_Eta", 100, -3, 3)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_Phi", "leadJet_Phi", 100, -3.14, 3.14)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_E", "leadJet_E", 200, 0, 1000)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_nConstituent", "leadJet_nConstituent", 200, 0, 200)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_JetRecHit_auxTDC3_frac0", "leadJet_JetRecHit_auxTDC3_frac0", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_JetRecHit_auxTDC3_frac1", "leadJet_JetRecHit_auxTDC3_frac1", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_JetRecHit_auxTDC3_frac2", "leadJet_JetRecHit_auxTDC3_frac2", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_JetRecHit_auxTDC3_frac3", "leadJet_JetRecHit_auxTDC3_frac3", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_JetRecHit_auxTDC3_energyWeightedFrac0", "leadJet_JetRecHit_auxTDC3_energyWeightedFrac0", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_JetRecHit_auxTDC3_energyWeightedFrac1", "leadJet_JetRecHit_auxTDC3_energyWeightedFrac1", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_JetRecHit_auxTDC3_energyWeightedFrac2", "leadJet_JetRecHit_auxTDC3_energyWeightedFrac2", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_JetRecHit_auxTDC3_energyWeightedFrac3", "leadJet_JetRecHit_auxTDC3_energyWeightedFrac3", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_JetRecHit_energy_avg", "leadJet_JetRecHit_energy_avg", 100, 0, 4)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_JetRecHit_energy_max", "leadJet_JetRecHit_energy_max", 400, 0, 200)
    hists = createHisto1D(hists, dfHB, "leadJet_HB_JetRecHit_energy_min", "leadJet_JetRecHit_energy_min", 100, 0, 2)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_Pt", "leadJet_Pt", 200, 0, 1000)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_Eta", "leadJet_Eta", 100, -3, 3)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_Phi", "leadJet_Phi", 100, -3.14, 3.14)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_E", "leadJet_E", 200, 0, 1000)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_nConstituent", "leadJet_nConstituent", 200, 0, 200)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_JetRecHit_auxTDC3_frac0", "leadJet_JetRecHit_auxTDC3_frac0", 100, 0, 1)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_JetRecHit_auxTDC3_frac1", "leadJet_JetRecHit_auxTDC3_frac1", 100, 0, 1)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_JetRecHit_auxTDC3_frac2", "leadJet_JetRecHit_auxTDC3_frac2", 100, 0, 1)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_JetRecHit_auxTDC3_frac3", "leadJet_JetRecHit_auxTDC3_frac3", 100, 0, 1)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_JetRecHit_auxTDC3_energyWeightedFrac0", "leadJet_JetRecHit_auxTDC3_energyWeightedFrac0", 100, 0, 1)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_JetRecHit_auxTDC3_energyWeightedFrac1", "leadJet_JetRecHit_auxTDC3_energyWeightedFrac1", 100, 0, 1)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_JetRecHit_auxTDC3_energyWeightedFrac2", "leadJet_JetRecHit_auxTDC3_energyWeightedFrac2", 100, 0, 1)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_JetRecHit_auxTDC3_energyWeightedFrac3", "leadJet_JetRecHit_auxTDC3_energyWeightedFrac3", 100, 0, 1)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_JetRecHit_energy_avg", "leadJet_JetRecHit_energy_avg", 100, 0, 4)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_JetRecHit_energy_max", "leadJet_JetRecHit_energy_max", 400, 0, 200)
    hists = createHisto1D(hists, dfHE, "leadJet_HE_JetRecHit_energy_min", "leadJet_JetRecHit_energy_min", 100, 0, 2)
    hists = createHisto1D(hists, df, "recHit_E_max", "RecHit_E_max", 400, 0, 200)
    hists = createHisto1D(hists, df, "recHit_E_min", "RecHit_E_min", 100, 0, 2)
    hists = createHisto1D(hists, df, "recHit_E_median", "RecHit_E_median", 100, 0, 4)
    hists = createHisto1D(hists, df, "recHit_E_mean", "RecHit_E_mean", 100, 0, 4)
    hists = createHisto1D(hists, df, "recHit_auxTDC3_frac0", "RecHit_auxTDC3_frac0", 200, -1, 1)
    hists = createHisto1D(hists, df, "recHit_auxTDC3_frac1", "RecHit_auxTDC3_frac1", 200, -1, 1)
    hists = createHisto1D(hists, df, "recHit_auxTDC3_frac2", "RecHit_auxTDC3_frac2", 200, -1, 1)
    hists = createHisto1D(hists, df, "recHit_auxTDC3_frac3", "RecHit_auxTDC3_frac3", 200, -1, 1)
    hists = createHisto1D(hists, df, "recHit_auxTDC3_energyWeightedFrac0", "RecHit_auxTDC3_energyWeightedFrac0", 200, -1, 1)
    hists = createHisto1D(hists, df, "recHit_auxTDC3_energyWeightedFrac1", "RecHit_auxTDC3_energyWeightedFrac1", 200, -1, 1)
    hists = createHisto1D(hists, df, "recHit_auxTDC3_energyWeightedFrac2", "RecHit_auxTDC3_energyWeightedFrac2", 200, -1, 1)
    hists = createHisto1D(hists, df, "recHit_auxTDC3_energyWeightedFrac3", "RecHit_auxTDC3_energyWeightedFrac3", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "recHit_HB_E_max", "RecHit_E_max", 400, 0, 200)
    hists = createHisto1D(hists, dfHB, "recHit_HB_E_min", "RecHit_E_min", 100, 0, 2)
    hists = createHisto1D(hists, dfHB, "recHit_HB_E_median", "RecHit_E_median", 100, 0, 4)
    hists = createHisto1D(hists, dfHB, "recHit_HB_E_mean", "RecHit_E_mean", 100, 0, 4)
    hists = createHisto1D(hists, dfHB, "recHit_HB_auxTDC3_frac0", "RecHit_auxTDC3_frac0", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "recHit_HB_auxTDC3_frac1", "RecHit_auxTDC3_frac1", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "recHit_HB_auxTDC3_frac2", "RecHit_auxTDC3_frac2", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "recHit_HB_auxTDC3_frac3", "RecHit_auxTDC3_frac3", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "recHit_HB_auxTDC3_energyWeightedFrac0", "RecHit_auxTDC3_energyWeightedFrac0", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "recHit_HB_auxTDC3_energyWeightedFrac1", "RecHit_auxTDC3_energyWeightedFrac1", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "recHit_HB_auxTDC3_energyWeightedFrac2", "RecHit_auxTDC3_energyWeightedFrac2", 200, -1, 1)
    hists = createHisto1D(hists, dfHB, "recHit_HB_auxTDC3_energyWeightedFrac3", "RecHit_auxTDC3_energyWeightedFrac3", 200, -1, 1)
    hists = createHisto1D(hists, dfHE, "recHit_HE_E_max", "RecHit_E_max", 400, 0, 200)
    hists = createHisto1D(hists, dfHE, "recHit_HE_E_min", "RecHit_E_min", 100, 0, 2)
    hists = createHisto1D(hists, dfHE, "recHit_HE_E_median", "RecHit_E_median", 100, 0, 4)
    hists = createHisto1D(hists, dfHE, "recHit_HE_E_mean", "RecHit_E_mean", 100, 0, 4)
    hists = createHisto1D(hists, dfHE, "recHit_HE_auxTDC3_frac0", "RecHit_auxTDC3_frac0", 200, -1, 1)
    hists = createHisto1D(hists, dfHE, "recHit_HE_auxTDC3_frac1", "RecHit_auxTDC3_frac1", 200, -1, 1)
    hists = createHisto1D(hists, dfHE, "recHit_HE_auxTDC3_frac2", "RecHit_auxTDC3_frac2", 200, -1, 1)
    hists = createHisto1D(hists, dfHE, "recHit_HE_auxTDC3_frac3", "RecHit_auxTDC3_frac3", 200, -1, 1)
    hists = createHisto1D(hists, dfHE, "recHit_HE_auxTDC3_energyWeightedFrac0", "RecHit_auxTDC3_energyWeightedFrac0", 200, -1, 1)
    hists = createHisto1D(hists, dfHE, "recHit_HE_auxTDC3_energyWeightedFrac1", "RecHit_auxTDC3_energyWeightedFrac1", 200, -1, 1)
    hists = createHisto1D(hists, dfHE, "recHit_HE_auxTDC3_energyWeightedFrac2", "RecHit_auxTDC3_energyWeightedFrac2", 200, -1, 1)
    hists = createHisto1D(hists, dfHE, "recHit_HE_auxTDC3_energyWeightedFrac3", "RecHit_auxTDC3_energyWeightedFrac3", 200, -1, 1)

    for metric in ["avg", "pTWeightedAvg", "med", "min", "max"]:
        for i in range(1,5):
            hists = createHisto2D(hists, dfHB, f"leadJet_HB_{metric}ConstituentHcalDepthEF{i}", f"leadJet_{metric}ConstituentHcalDepthEF{i}", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
        for i in range(1,8):
            hists = createHisto2D(hists, dfHE, f"leadJet_HE_{metric}ConstituentHcalDepthEF{i}", f"leadJet_{metric}ConstituentHcalDepthEF{i}", 50, 0, 1, "leadJet_Pt", 20, 0, 1000)
        hists = createHisto2D(hists, dfHB, f"leadJet_HB_{metric}ConstituentRawHoverE", f"leadJet_{metric}ConstituentRawHoverE", 200, 0, 100, "leadJet_Pt", 20, 0, 1000)
        hists = createHisto2D(hists, dfHB, f"leadJet_HB_{metric}ConstituentHoverE", f"leadJet_{metric}ConstituentHoverE", 200, 0, 100, "leadJet_Pt", 20, 0, 1000)
        hists = createHisto2D(hists, dfHE, f"leadJet_HE_{metric}ConstituentRawHoverE", f"leadJet_{metric}ConstituentRawHoverE", 200, 0, 100, "leadJet_Pt", 20, 0, 1000)
        hists = createHisto2D(hists, dfHE, f"leadJet_HE_{metric}ConstituentHoverE", f"leadJet_{metric}ConstituentHoverE", 200, 0, 100, "leadJet_Pt", 20, 0, 1000)
    
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
