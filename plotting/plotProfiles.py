###################################################################################################
#   plotProfiles.py                                                                               #
#   Description: create profile plots of jet constituent info from histograms created by          #
#   writeHistograms.py                                                                            #
#   Author: Elliott Kauffman                                                                      #
###################################################################################################

###################################################################################################
# IMPORTS

import argparse
import datetime
import matplotlib.pyplot as plt
import mplhep as hep
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import uproot

###################################################################################################
# GLOBAL VARIABLES

def constructSampleDict():

    mMed_values = [100, 250, 500, 750, 1000, 1500]
    mDark_values = [10, 20]
    ctau_values = [1, 100, 1000, 1500, 2000]
    
    sampleDict = {}
    
    for mMed in mMed_values:
        for mDark in mDark_values:
            for ctau in ctau_values:
                key = f"EMJ_s-channel_mMed-{mMed}_mDark-{mDark}_ctau-{ctau}_unflavored-down"
                sampleDict[key] = {
                    "histFile": f"histograms/histograms_EMJ_s-channel_mMed-{mMed}_mDark-{mDark}_ctau-{ctau}_unflavored-down.root",
                    "printName": f"EMJ (mMed={mMed}, mDark={mDark}, ctau={ctau})"
                }
                
    sampleDict["QCD_Combined"] = {
        "histFile": "histograms/histograms_QCD_Combined.root",
        "printName": "QCD",
    }
    
    qcd_bins = [
        "30to50", "50to80", "80to120", "120to170", "170to300", "300to470", "470to600", "600to800", "800to1000", "1000"
    ]
    
    for bin in qcd_bins:
        key = f"QCD_Bin-Pt-{bin}"
        if bin == "1000":
            print_name = "QCD (pT>1000 GeV)"
        else:
            pt_range = bin.replace("to", " to ")
            print_name = f"QCD (pT={pt_range} GeV)"
        
        sampleDict[key] = {
            "histFile": f"histograms/histograms_QCD_Bin-Pt-{bin}.root",
            "printName": print_name
        }

    print(sampleDict)
    
    return sampleDict


###################################################################################################
# METHODS

def parseArgs() -> argparse.Namespace:

    p = argparse.ArgumentParser()
    p.add_argument(
        "--outdir",
        "-o",
        help="Name of the directory where plots are stored.",
        default="./plots/",
    )
    return p.parse_args()
    
def plotHCalDepthEFProfiles(sampleDict, histograms, samples, y_label, outdir, plot_filename, region, pt_min=None, pt_max=None):

    colors = ['black', 'crimson', 'slateblue', 'tomato', 'seagreen', 'teal', 'darkmagenta']

    if region not in ["HE", "HB"]: return
    
    if region=="HE":
        hist_names = ["1", "2", "3", "4", "5", "6", "7"]
        depth_edges = list(range(8))
        plotting_centers = np.arange(1, 8)-0.5
    else:
        hist_names = ["1", "2", "3", "4"]
        depth_edges = list(range(5))
        plotting_centers = np.arange(1, 5)-0.5
    
    bin_edges = histograms[0][0].axes[0].edges()
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    pt_bin_edges = histograms[0][0].axes[1].edges()
    if pt_min and pt_max:
        pt_start = np.searchsorted(pt_bin_edges, pt_min, side="left") - 1
        pt_end = np.searchsorted(pt_bin_edges, pt_max, side="right") - 1
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    
    for i in range(len(histograms)):
        if pt_min and pt_max:
            values_2d = np.array(
                [np.sum(hist.values()[:, pt_start:pt_end], axis=1) for hist in histograms[i]]
            ).T
        else:
            values_2d = np.array(
                [np.sum(hist.values(), axis=1) for hist in histograms[i]]
            ).T
            
        values_2d_for_profile = values_2d[1:]  # exclude 0-th bin
        bin_centers_for_profile = bin_centers[1:]
        sum_weights = np.sum(values_2d_for_profile, axis=0)
        sum_weighted_bins = np.sum(values_2d_for_profile * bin_centers_for_profile[:, None], axis=0)
        profile = np.divide(sum_weighted_bins, sum_weights, where=(sum_weights != 0))
        profile = np.nan_to_num(profile, nan=0.0, posinf=0.0, neginf=0.0)
        profile[(profile < 0) | (profile > 1)] = 0

        sum_weighted_bins_sq = np.sum(values_2d_for_profile * (bin_centers_for_profile[:, None])**2, axis=0)
        variance = np.divide(sum_weighted_bins_sq, sum_weights, where=(sum_weights != 0)) - profile**2
        variance = np.maximum(variance, 0)  # Ensure variance is non-negative
        stddev = np.sqrt(variance)
        stderr = np.divide(stddev, np.sqrt(sum_weights), where=(sum_weights != 0))
        stderr = np.nan_to_num(stderr, nan=0.0, posinf=0.0, neginf=0.0)

        profile_errors_down = np.maximum(0, np.minimum(stderr, profile))  # Prevent going below 0
        profile_errors_up = np.maximum(0, np.minimum(stderr, np.where(profile == 0, 0, 1 - profile)))  # Prevent exceeding 1
        asymmetric_errors = [profile_errors_down, profile_errors_up]

        
        ax.errorbar(plotting_centers,
                    profile,
                    yerr=asymmetric_errors,
                    fmt='-o',
                    color=colors[i],
                    capsize=5
                    )
        
    ax.set_xticks(np.arange(len(hist_names)) + 0.5)
    ax.legend([sampleDict[samples[i]]["printName"] for i in range(len(samples))], prop={'size': 18}, loc="upper left")
    ax.set_xlabel("HCAL Depth Level")
    ax.set_ylim([0,1])
    ax.set_ylabel(y_label, fontsize=18)
    ax.set_xticklabels(hist_names)
    ax.grid(True, linestyle="--", alpha=0.5)
    
    fig.savefig(outdir + plot_filename + ".png")
    plt.close(fig)
        
    
###################################################################################################
# MAIN METHOD

def main():

    args = parseArgs()
    todaysDate = datetime.date.today().strftime('%Y%m%d')
    hep.style.use("CMS")
    
    sampleDict = constructSampleDict()
    
    pt_min_vals = [0,50,150,300]
    pt_max_vals = [50,150,300,500]
    categories = ['avg', 'med', 'min', 'max','pTWeightedAvg']
    category_names = ['Average', 'Median', 'Minimum', 'Maximum', 'pT-Weighted Average']
    regions = ['HE', 'HB']
    mDark_vals = [10,20]
    ctau_vals = [1,100,1000,1500,2000]
    
    for k, category in enumerate(categories):
        print("category = ", category)
        
        for mDark_val in mDark_vals:
            print("    mDark = ", mDark_val)
            
            for ctau_val in ctau_vals:
                print("        ctau = ", ctau_val)
                samples = [
                    "QCD_Combined",
                    f"EMJ_s-channel_mMed-100_mDark-{mDark_val}_ctau-{ctau_val}_unflavored-down",
                    f"EMJ_s-channel_mMed-250_mDark-{mDark_val}_ctau-{ctau_val}_unflavored-down",
                    f"EMJ_s-channel_mMed-500_mDark-{mDark_val}_ctau-{ctau_val}_unflavored-down",
                    f"EMJ_s-channel_mMed-750_mDark-{mDark_val}_ctau-{ctau_val}_unflavored-down",
                    f"EMJ_s-channel_mMed-1000_mDark-{mDark_val}_ctau-{ctau_val}_unflavored-down",
                    f"EMJ_s-channel_mMed-1500_mDark-{mDark_val}_ctau-{ctau_val}_unflavored-down",
                ]

                files = [uproot.open(sampleDict[samples[i]]["histFile"]) for i in range(len(samples))]
                
                for region in regions:
                    print("            region = ", region)
                    
                    if region=="HE": range_end = 8
                    else: range_end = 5
                    
                    histograms = [
                        [
                            files[j][f"leadJet_{region}_{category}ConstituentHcalDepthEF{i}"] for i in range(1,range_end)
                        ] for j in range(len(samples))
                    ]
                    
                    plotHCalDepthEFProfiles(
                        sampleDict,
                        histograms,
                        samples,
                        f"Profile: Lead Jet {category_names[k]} HCAL Energy Fraction",
                        args.outdir,
                        f"leadJet_{region}_{category}ConstituentHcalDepthEF_mDark{mDark_val}_ctau{ctau_val}_allPt_{todaysDate}",
                        region
                    )
                    for i in range(len(pt_min_vals)):
                        plotHCalDepthEFProfiles(
                            sampleDict,
                            histograms,
                            samples,
                            f"Profile: Lead Jet {category_names[k]} HCAL Energy Fraction",
                            args.outdir,
                            f"leadJet_{region}_{category}ConstituentHcalDepthEF_mDark{mDark_val}_ctau{ctau_val}_Pt{pt_min_vals[i]}to{pt_max_vals[i]}_{todaysDate}",
                            region,
                            pt_min=pt_min_vals[i],
                            pt_max=pt_max_vals[i],
                        )
                        
    mMed_vals = [100,250,500,750,1000,1500]
    
    for k, category in enumerate(categories):
        print("category = ", category)
        
        for mDark_val in mDark_vals:
            print("    mDark = ", mDark_val)
            
            for mMed_val in mMed_vals:
                print("        mMed = ", mMed_val)
                samples = [
                    "QCD_Combined",
                    f"EMJ_s-channel_mMed-{mMed_val}_mDark-{mDark_val}_ctau-1_unflavored-down",
                    f"EMJ_s-channel_mMed-{mMed_val}_mDark-{mDark_val}_ctau-100_unflavored-down",
                    f"EMJ_s-channel_mMed-{mMed_val}_mDark-{mDark_val}_ctau-1000_unflavored-down",
                    f"EMJ_s-channel_mMed-{mMed_val}_mDark-{mDark_val}_ctau-1500_unflavored-down",
                    f"EMJ_s-channel_mMed-{mMed_val}_mDark-{mDark_val}_ctau-2000_unflavored-down",
                ]
                
                
                files = [uproot.open(sampleDict[samples[i]]["histFile"]) for i in range(len(samples))]
                
                for region in regions:
                    print("            region = ", region)
                        
                    if region=="HE": range_end = 8
                    else: range_end = 5
                    
                    histograms = [
                        [
                            files[j][f"leadJet_{region}_{category}ConstituentHcalDepthEF{i}"] for i in range(1,range_end)
                        ] for j in range(len(samples))
                    ]
                    
                    plotHCalDepthEFProfiles(
                        sampleDict,
                        histograms,
                        samples,
                        f"Profile: Lead Jet {category_names[k]} HCAL Energy Fraction",
                        args.outdir,
                        f"leadJet_{region}_{category}ConstituentHcalDepthEF_mMed{mMed_val}_mDark{mDark_val}_allPt_{todaysDate}",
                        region
                    )
                    for i in range(len(pt_min_vals)):
                        plotHCalDepthEFProfiles(
                            sampleDict,
                            histograms,
                            samples,
                            f"Profile: Lead Jet {category_names[k]} HCAL Energy Fraction",
                            args.outdir,
                            f"leadJet_{region}_{category}ConstituentHcalDepthEF_mMed{mMed_val}_mDark{mDark_val}_Pt{pt_min_vals[i]}to{pt_max_vals[i]}_{todaysDate}",
                            region,
                            pt_min=pt_min_vals[i],
                            pt_max=pt_max_vals[i],
                        )



###################################################################################################
# RUN SCRIPT

if __name__ == "__main__":
    main()
