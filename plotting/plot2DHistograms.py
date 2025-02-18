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
# METHODS

def parseArgs() -> argparse.Namespace:

    p = argparse.ArgumentParser()
    p.add_argument(
        "--input",
        "-i",
        help="Name of the ROOT file containing the histograms to plot.",
    )
    p.add_argument(
        "--dataset",
        "-d",
        help="Name of the dataset we are creating plots for",
    )
    p.add_argument(
        "--outdir",
        "-o",
        help="Name of the directory where plots are stored.",
        default="./plots/",
    )
    return p.parse_args()
    
def plotHCalDepthEFs(histograms, y_label, outdir, plot_filename, norm_method='log'):

    hist_names = ["1", "2", "3", "4", "5", "6", "7"]
    depth_edges = list(range(8))
    bin_edges = histograms[0].axes[0].edges()
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    X, Y = np.meshgrid(depth_edges, bin_edges)
    values_2d = np.array([hist.values() for hist in histograms]).T
    masked_vals = np.ma.masked_where(values_2d == 0, values_2d)
    
    sum_weights = np.sum(values_2d, axis=0)  # Sum of counts per depth
    sum_weighted_bins = np.sum(values_2d * bin_centers[:, None], axis=0)  # Weighted sum per depth
    profile = np.divide(sum_weighted_bins, sum_weights, where=(sum_weights != 0))
    
    variance_numerator = np.sum(values_2d * (bin_centers[:, None] - profile) ** 2, axis=0)
    weighted_std_dev = np.sqrt(np.divide(variance_numerator, sum_weights, where=(sum_weights != 0)))
    
    lower_err = np.maximum(profile - weighted_std_dev, 0)
    upper_err = np.minimum(profile + weighted_std_dev, 1)
    asymmetric_errors = [profile - lower_err, upper_err - profile]
    
    fig, axs = plt.subplots(2, 1, figsize=(10, 12),
                            gridspec_kw={'height_ratios': [3, 1]},
                            sharex=True)
                            
    ax_pcolormesh = axs[0]  # Upper plot
    ax_profile = axs[1]     # Lower profile plot
    
    fig.subplots_adjust(right=0.8)
    cax = fig.add_axes([0.85, 0.36, 0.05, 0.51])
    
    pcm = ax_pcolormesh.pcolormesh(X, Y, masked_vals, shading="auto", cmap="plasma", norm=norm_method)
    ax_pcolormesh.set_ylabel(y_label, fontsize=18)
    ax_pcolormesh.set_xticks(np.arange(len(hist_names)) + 0.5)
    
    fig.colorbar(pcm, cax=cax, orientation='vertical')
    
    ax_profile.errorbar(np.arange(1, 8)-0.5, profile, yerr=asymmetric_errors, fmt='o', color='#5790fc', capsize=5)
    ax_profile.set_xlabel("HCAL Depth Level")
    ax_profile.set_ylabel("Profile", fontsize=18)
    ax_profile.set_xticklabels(hist_names)
    ax_profile.grid(True, linestyle="--", alpha=0.5)
    
    #fig.tight_layout()
    #fig.subplots_adjust(wspace=0.05, hspace=0.05)
    fig.savefig(outdir + plot_filename + ".pdf")
    fig.savefig(outdir + plot_filename + ".png")
    
###################################################################################################
# MAIN METHOD

def main():

    args = parseArgs()
    
    todaysDate = datetime.date.today().strftime('%Y%m%d')
    
    f = uproot.open(args.input)
    
    hep.style.use("CMS")
    
    # create leading jet average hcal energy fraction plot
    histograms = [f[f"leadJet_avgConstituentHcalDepthEF{i}"] for i in range(1,8)]
    plotHCalDepthEFs(histograms,
                     "Lead Jet Average HCAL Energy Fraction",
                     args.outdir,
                     f"leadJet_avgConstituentHcalDepthEF_{args.dataset}_{todaysDate}")
                     
    # create leading jet pt-weighted average hcal energy fraction plot
    histograms = [f[f"leadJet_pTWeightedAvgConstituentHcalDepthEF{i}"] for i in range(1,8)]
    plotHCalDepthEFs(histograms,
                     "Lead Jet pT-Weighted Average HCAL Energy Fraction",
                     args.outdir,
                     f"leadJet_pTWeightedAvgConstituentHcalDepthEF_{args.dataset}_{todaysDate}")
                     
    # create leading jet median hcal energy fraction plot
    histograms = [f[f"leadJet_medConstituentHcalDepthEF{i}"] for i in range(1,8)]
    plotHCalDepthEFs(histograms,
                     "Lead Jet Median HCAL Energy Fraction",
                     args.outdir,
                     f"leadJet_medConstituentHcalDepthEF_{args.dataset}_{todaysDate}")
                     
    # create leading jet minimum hcal energy fraction plot
    histograms = [f[f"leadJet_minConstituentHcalDepthEF{i}"] for i in range(1,8)]
    plotHCalDepthEFs(histograms,
                     "Lead Jet Minimum HCAL Energy Fraction",
                     args.outdir,
                     f"leadJet_minConstituentHcalDepthEF_{args.dataset}_{todaysDate}")
                     
    # create leading jet maximum hcal energy fraction plot
    histograms = [f[f"leadJet_maxConstituentHcalDepthEF{i}"] for i in range(1,8)]
    plotHCalDepthEFs(histograms,
                     "Lead Jet Maximum HCAL Energy Fraction",
                     args.outdir,
                     f"leadJet_maxConstituentHcalDepthEF_{args.dataset}_{todaysDate}")
    

###################################################################################################
# RUN SCRIPT

if __name__ == "__main__":
    main()
