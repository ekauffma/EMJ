###################################################################################################
#   plot2DHistograms.py                                                                           #
#   Description: create 2d plots of jet constituent info from histograms created by               #
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

inputFileDict = {
    "QCD_Combined": {
        "histFile": "histograms/histograms_QCD_Combined.root",
        "printName": "QCD",
        "color": "black"
    },
    "QCD_Bin-Pt-30to50": {
        "histFile": "histograms/histograms_QCD_Bin-Pt-30to50.root",
        "printName": "QCD (pT=30 to 50 GeV)",
        "color": "lightcoral"
    },
    "QCD_Bin-Pt-50to80": {
        "histFile": "histograms/histograms_QCD_Bin-Pt-50to80.root",
        "printName": "QCD (pT=50 to 80 GeV)",
        "color": "tomato"
    },
    "QCD_Bin-Pt-80to120": {
        "histFile": "histograms/histograms_QCD_Bin-Pt-80to120.root",
        "printName": "QCD (pT=80 to 120 GeV)",
        "color": "chocolate"
    },
    "QCD_Bin-Pt-120to170": {
        "histFile": "histograms/histograms_QCD_Bin-Pt-120to170.root",
        "printName": "QCD (pT=120 to 170 GeV)",
        "color": "orange"
    },
    "QCD_Bin-Pt-170to300": {
        "histFile": "histograms/histograms_QCD_Bin-Pt-170to300.root",
        "printName": "QCD (pT=170 to 300 GeV)",
        "color": "gold"
    },
    "QCD_Bin-Pt-300to470": {
        "histFile": "histograms/histograms_QCD_Bin-Pt-300to470.root",
        "printName": "QCD (pT=300 to 470 GeV)",
        "color": "olive"
    },
    "QCD_Bin-Pt-470to600": {
        "histFile": "histograms/histograms_QCD_Bin-Pt-470to600.root",
        "printName": "QCD (pT=470 to 600 GeV)",
        "color": "lime"
    },
    "QCD_Bin-Pt-600to800": {
        "histFile": "histograms/histograms_QCD_Bin-Pt-600to800.root",
        "printName": "QCD (pT=600 to 800 GeV)",
        "color": "darkturquoise"
    },
    "QCD_Bin-Pt-800to1000": {
        "histFile": "histograms/histograms_QCD_Bin-Pt-800to1000.root",
        "printName": "QCD (pT=800 to 1000 GeV)",
        "color": "royalblue"
    },
    "QCD_Bin-Pt-1000": {
        "histFile": "histograms/histograms_QCD_Bin-Pt-1000.root",
        "printName": "QCD (pT>1000 GeV)",
        "color": "violet"
    },
    "EMJ_s-channel_mMed-100_mDark-10_ctau-1000_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-10_ctau-1000_unflavored-down.root",
        "printName": "EMJ (mMed=100, mDark=10, ctau=1000)",
        "color": "indianred"
    },
    "EMJ_s-channel_mMed-100_mDark-10_ctau-100_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-10_ctau-100_unflavored-down.root",
        "printName": "EMJ (mMed=100, mDark=10, ctau=100)",
        "color": "saddlebrown"
    },
    "EMJ_s-channel_mMed-100_mDark-10_ctau-1_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-10_ctau-1_unflavored-down.root",
        "printName": "EMJ (mMed=100, mDark=10, ctau=1)",
        "color": "darkgoldenrod"
    },
    "EMJ_s-channel_mMed-100_mDark-20_ctau-1000_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-20_ctau-1000_unflavored-down.root",
        "printName": "EMJ (mMed=100, mDark=20, ctau=1000)",
        "color": "seagreen"
    },
    "EMJ_s-channel_mMed-100_mDark-10_ctau-100_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-10_ctau-100_unflavored-down.root",
        "printName": "EMJ (mMed=100, mDark=20, ctau=100)",
        "color": "steelblue"
    },
    "EMJ_s-channel_mMed-100_mDark-20_ctau-1_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-20_ctau-1_unflavored-down.root",
        "printName": "EMJ (mMed=100, mDark=20, ctau=1)",
        "color": "mediumorchid"
    },
    "EMJ_s-channel_mMed-250_mDark-10_ctau-1000_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-10_ctau-1000_unflavored-down.root",
        "printName": "EMJ (mMed=250, mDark=10, ctau=1000)",
        "color": "firebrick"
    },
    "EMJ_s-channel_mMed-250_mDark-10_ctau-100_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-10_ctau-100_unflavored-down.root",
        "printName": "EMJ (mMed=250, mDark=10, ctau=100)",
        "color": "darkorange"
    },
    "EMJ_s-channel_mMed-250_mDark-10_ctau-1_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-10_ctau-1_unflavored-down.root",
        "printName": "EMJ (mMed=250, mDark=10, ctau=1)",
        "color": "yellowgreen"
    },
    "EMJ_s-channel_mMed-250_mDark-20_ctau-1000_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-20_ctau-1000_unflavored-down.root",
        "printName": "EMJ (mMed=250, mDark=20, ctau=1000)",
        "color": "darkcyan"
    },
    "EMJ_s-channel_mMed-250_mDark-10_ctau-100_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-10_ctau-100_unflavored-down.root",
        "printName": "EMJ (mMed=250, mDark=20, ctau=100)",
        "color": "blue"
    },
    "EMJ_s-channel_mMed-250_mDark-20_ctau-1_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-20_ctau-1_unflavored-down.root",
        "printName": "EMJ (mMed=250, mDark=20, ctau=1)",
        "color": "indigo"
    },
    "EMJ_s-channel_mMed-500_mDark-10_ctau-1000_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-500_mDark-10_ctau-1000_unflavored-down.root",
        "printName": "EMJ (mMed=500, mDark=10, ctau=1000)",
        "color": "red"
    },
    "EMJ_s-channel_mMed-500_mDark-10_ctau-100_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-500_mDark-10_ctau-100_unflavored-down.root",
        "printName": "EMJ (mMed=500, mDark=10, ctau=100)",
        "color": "lightsalmon"
    },
    "EMJ_s-channel_mMed-500_mDark-10_ctau-1_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-500_mDark-10_ctau-1_unflavored-down.root",
        "printName": "EMJ (mMed=500, mDark=10, ctau=1)",
        "color": "forestgreen"
    },
    "EMJ_s-channel_mMed-500_mDark-20_ctau-1000_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-500_mDark-20_ctau-1000_unflavored-down.root",
        "printName": "EMJ (mMed=500, mDark=20, ctau=1000)",
        "color": "mediumturquoise"
    },
    "EMJ_s-channel_mMed-500_mDark-10_ctau-100_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-500_mDark-10_ctau-100_unflavored-down.root",
        "printName": "EMJ (mMed=500, mDark=20, ctau=100)",
        "color": "mediumblue"
    },
    "EMJ_s-channel_mMed-500_mDark-20_ctau-1_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-500_mDark-20_ctau-1_unflavored-down.root",
        "printName": "EMJ (mMed=500, mDark=20, ctau=1)",
        "color": "purple"
    },
    "EMJ_s-channel_mMed-750_mDark-10_ctau-1000_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-750_mDark-10_ctau-1000_unflavored-down.root",
        "printName": "EMJ (mMed=750, mDark=10, ctau=1000)",
        "color": "orangered"
    },
    "EMJ_s-channel_mMed-750_mDark-10_ctau-100_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-750_mDark-10_ctau-100_unflavored-down.root",
        "printName": "EMJ (mMed=750, mDark=10, ctau=100)",
        "color": "goldenrod"
    },
    "EMJ_s-channel_mMed-750_mDark-10_ctau-1_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-750_mDark-10_ctau-1_unflavored-down.root",
        "printName": "EMJ (mMed=750, mDark=10, ctau=1)",
        "color": "olivedrab"
    },
    "EMJ_s-channel_mMed-750_mDark-20_ctau-1000_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-750_mDark-20_ctau-1000_unflavored-down.root",
        "printName": "EMJ (mMed=750, mDark=20, ctau=1000)",
        "color": "aquamarine"
    },
    "EMJ_s-channel_mMed-750_mDark-10_ctau-100_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-750_mDark-10_ctau-100_unflavored-down.root",
        "printName": "EMJ (mMed=750, mDark=20, ctau=100)",
        "color": "rebeccapurple"
    },
    "EMJ_s-channel_mMed-750_mDark-20_ctau-1_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-750_mDark-20_ctau-1_unflavored-down.root",
        "printName": "EMJ (mMed=750, mDark=20, ctau=1)",
        "color": "deeppink"
    },
    "EMJ_s-channel_mMed-1000_mDark-10_ctau-1000_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-1000_mDark-10_ctau-1000_unflavored-down.root",
        "printName": "EMJ (mMed=1000, mDark=10, ctau=1000)",
        "color": "coral"
    },
    "EMJ_s-channel_mMed-1000_mDark-10_ctau-100_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-1000_mDark-10_ctau-100_unflavored-down.root",
        "printName": "EMJ (mMed=1000, mDark=10, ctau=100)",
        "color": "orange"
    },
    "EMJ_s-channel_mMed-1000_mDark-10_ctau-1_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-1000_mDark-10_ctau-1_unflavored-down.root",
        "printName": "EMJ (mMed=1000, mDark=10, ctau=1)",
        "color": "green"
    },
    "EMJ_s-channel_mMed-1000_mDark-20_ctau-1000_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-1000_mDark-20_ctau-1000_unflavored-down.root",
        "printName": "EMJ (mMed=1000, mDark=20, ctau=1000)",
        "color": "teal"
    },
    "EMJ_s-channel_mMed-1000_mDark-10_ctau-100_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-1000_mDark-10_ctau-100_unflavored-down.root",
        "printName": "EMJ (mMed=1000, mDark=20, ctau=100)",
        "color": "navy"
    },
    "EMJ_s-channel_mMed-1000_mDark-20_ctau-1_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-1000_mDark-20_ctau-1_unflavored-down.root",
        "printName": "EMJ (mMed=1000, mDark=20, ctau=1)",
        "color": "darkviolet"
    },
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-1000_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-1500_mDark-10_ctau-1000_unflavored-down.root",
        "printName": "EMJ (mMed=1500, mDark=10, ctau=1000)",
        "color": "salmon"
    },
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-100_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-1500_mDark-10_ctau-100_unflavored-down.root",
        "printName": "EMJ (mMed=1500, mDark=10, ctau=100)",
        "color": "peru"
    },
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-1_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-1500_mDark-10_ctau-1_unflavored-down.root",
        "printName": "EMJ (mMed=1500, mDark=10, ctau=1)",
        "color": "mediumseagreen"
    },
    "EMJ_s-channel_mMed-1500_mDark-20_ctau-1000_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-1500_mDark-20_ctau-1000_unflavored-down.root",
        "printName": "EMJ (mMed=1500, mDark=20, ctau=1000)",
        "color": "deepskyblue"
    },
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-100_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-1500_mDark-10_ctau-100_unflavored-down.root",
        "printName": "EMJ (mMed=1500, mDark=20, ctau=100)",
        "color": "slateblue"
    },
    "EMJ_s-channel_mMed-1500_mDark-20_ctau-1_unflavored-down": {
        "histFile": "histograms/histograms_EMJ_s-channel_mMed-1500_mDark-20_ctau-1_unflavored-down.root",
        "printName": "EMJ (mMed=1500, mDark=20, ctau=1)",
        "color": "orchid"
    }
}

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
    
def plotHCalDepthEFProfiles(histograms, samples, y_label, outdir, plot_filename, pt_min=0, pt_max=1000):

    hist_names = ["1", "2", "3", "4", "5", "6", "7"]
    depth_edges = list(range(8))
    bin_edges = histograms[0][0].axes[0].edges()[1:] # skip the first bin
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    pt_bin_edges = histograms[0][0].axes[1].edges()
    pt_start = np.searchsorted(pt_bin_edges, pt_min, side="left") - 1
    pt_end = np.searchsorted(pt_bin_edges, pt_max, side="right") - 1
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    
    for i in range(len(histograms)):
        values_2d = np.array(
            [np.sum(hist.values()[:, pt_start:pt_end], axis=1)[1:] for hist in histograms[i]]
        ).T
        sum_weights = np.sum(values_2d, axis=0)
        sum_weighted_bins = np.sum(values_2d * bin_centers[:, None], axis=0)
        profile = np.divide(sum_weighted_bins, sum_weights, where=(sum_weights != 0))
        profile = np.nan_to_num(profile, nan=0.0, posinf=0.0, neginf=0.0)
        profile[(profile < 0) | (profile > 1)] = 0
        variance_numerator = np.sum(values_2d * (bin_centers[:, None] - profile) ** 2, axis=0)
        weighted_std_dev = np.sqrt(np.divide(variance_numerator, sum_weights, where=(sum_weights != 0)))
        weighted_std_dev = np.nan_to_num(weighted_std_dev, nan=0.0, posinf=0.0, neginf=0.0)
        profile_errors_down = np.minimum(weighted_std_dev, profile)  # Prevent going below 0
        profile_errors_up = np.maximum(0, np.minimum(weighted_std_dev, 1 - profile))  # Prevent exceeding 1
        asymmetric_errors = [profile_errors_down, profile_errors_up]
        
        ax.errorbar(np.arange(1, 8)-0.5,
                    profile,
                    yerr=asymmetric_errors,
                    fmt='o',
                    color=inputFileDict[samples[i]]["color"],
                    capsize=5
                    )
        
    ax.set_xticks(np.arange(len(hist_names)) + 0.5)
    ax.legend([inputFileDict[samples[i]]["printName"] for i in range(len(samples))], prop={'size': 18})
    ax.set_xlabel("HCAL Depth Level")
    ax.set_ylabel(y_label, fontsize=18)
    ax.set_xticklabels(hist_names)
    ax.grid(True, linestyle="--", alpha=0.5)
    
    fig.savefig(outdir + plot_filename + ".pdf")
    fig.savefig(outdir + plot_filename + ".png")
        
    
###################################################################################################
# MAIN METHOD

def main():

    args = parseArgs()
    todaysDate = datetime.date.today().strftime('%Y%m%d')
    hep.style.use("CMS")
    
    pt_min_vals = [0,50,100,150,200,250]
    pt_max_vals = [50,100,150,200,250,300]
    categories = ['avg', 'med', 'min', 'max','pTWeightedAvg']
    category_names = ['Average', 'Median', 'Minimum', 'Maximum', 'pT-Weighted Average']
    mDark_vals = [10,20]
    ctau_vals = [1,100,1000]
    
    for k, category in enumerate(categories):
        for mDark_val in mDark_vals:
            for ctau_val in ctau_vals:
                samples = [
                    "QCD_Combined",
                    f"EMJ_s-channel_mMed-100_mDark-{mDark_val}_ctau-{ctau_val}_unflavored-down",
                    f"EMJ_s-channel_mMed-250_mDark-{mDark_val}_ctau-{ctau_val}_unflavored-down",
                    f"EMJ_s-channel_mMed-500_mDark-{mDark_val}_ctau-{ctau_val}_unflavored-down",
                    f"EMJ_s-channel_mMed-750_mDark-{mDark_val}_ctau-{ctau_val}_unflavored-down",
                    f"EMJ_s-channel_mMed-1000_mDark-{mDark_val}_ctau-{ctau_val}_unflavored-down",
                    f"EMJ_s-channel_mMed-1500_mDark-{mDark_val}_ctau-{ctau_val}_unflavored-down",
                ]
                files = [uproot.open(inputFileDict[samples[i]]["histFile"]) for i in range(len(samples))]
                histograms = [
                    [
                        files[j][f"leadJet_{category}ConstituentHcalDepthEF{i}"] for i in range(1,8)
                    ] for j in range(len(samples))
                ]
                
                plotHCalDepthEFProfiles(
                    histograms,
                    samples,
                    f"Profile: Lead Jet {category_names[k]} HCAL Energy Fraction",
                    args.outdir,
                    f"leadJet_{category}ConstituentHcalDepthEF_mDark{mDark_val}_ctau{ctau_val}_allPt_{todaysDate}"
                )
                for i in range(len(pt_min_vals)):
                    plotHCalDepthEFProfiles(
                        histograms,
                        samples,
                        f"Profile: Lead Jet {category_names[k]} HCAL Energy Fraction",
                        args.outdir,
                        f"leadJet_{category}ConstituentHcalDepthEF_mDark{mDark_val}_ctau{ctau_val}_Pt{pt_min_vals[i]}to{pt_max_vals[i]}_{todaysDate}"
                        pt_min=pt_min_vals[i],
                        pt_max=pt_max_vals[i],
                    )

###################################################################################################
# RUN SCRIPT

if __name__ == "__main__":
    main()
