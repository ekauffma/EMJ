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
    
def plotHCalDepthEFs(histograms, y_label, outdir, plot_filename, region, norm_method='log', pt_min=None, pt_max=None):

    if region not in ["HE", "HB"]: return
    
    if region=="HE":
        hist_names = ["1", "2", "3", "4", "5", "6", "7"]
        depth_edges = np.array([1, 2, 5, 8, 11, 15, 19, 23])
        #depth_edges = list(range(8))
        plotting_centers = (depth_edges[:-1] + depth_edges[1:]) / 2
        #plotting_centers = np.arange(1, 8)-0.5
    else:
        hist_names = ["1", "2", "3", "4"]
        depth_edges = np.array([1, 2, 6, 11, 18])
        #depth_edges = list(range(5))
        plotting_centers = (depth_edges[:-1] + depth_edges[1:]) / 2
        #plotting_centers = np.arange(1, 5)-0.5
        
    bin_edges = histograms[0].axes[0].edges()
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    if pt_min and pt_max:
        pt_bin_edges = histograms[0].axes[1].edges()
        pt_start = np.searchsorted(pt_bin_edges, pt_min, side="left") - 1
        pt_end = np.searchsorted(pt_bin_edges, pt_max, side="right") - 1
    
    #X, Y = np.meshgrid(depth_edges[:-1], bin_edges)
    X, Y = np.meshgrid(depth_edges, bin_edges)
    if pt_min and pt_max:
        values_2d = np.array(
            [np.sum(hist.values()[:, pt_start:pt_end], axis=1) for hist in histograms]
        ).T
    else:
        values_2d = np.array(
            [np.sum(hist.values(), axis=1) for hist in histograms]
        ).T
    if sum(sum(values_2d))==0: return
    masked_vals = np.ma.masked_where(values_2d == 0, values_2d)
    
    sum_weights = np.sum(values_2d, axis=0)
    sum_weighted_bins = np.sum(values_2d * bin_centers[:, None], axis=0)
    profile = np.divide(sum_weighted_bins, sum_weights, where=(sum_weights != 0))
    profile = np.nan_to_num(profile, nan=0.0, posinf=0.0, neginf=0.0)
    profile[(profile < 0) | (profile > 1)] = 0

    sum_weighted_bins_sq = np.sum(values_2d * (bin_centers[:, None])**2, axis=0)
    variance = np.divide(sum_weighted_bins_sq, sum_weights, where=(sum_weights != 0)) - profile**2
    variance = np.maximum(variance, 0)  # Ensure variance is non-negative
    stddev = np.sqrt(variance)
    stderr = np.divide(stddev, np.sqrt(sum_weights), where=(sum_weights != 0))
    stderr = np.nan_to_num(stderr, nan=0.0, posinf=0.0, neginf=0.0)

    profile_errors_down = np.maximum(0, np.minimum(stderr, profile))  # Prevent going below 0
    profile_errors_up = np.maximum(0, np.minimum(stderr, 1 - profile))  # Prevent exceeding 1
    asymmetric_errors = [profile_errors_down, profile_errors_up]

    fig, axs = plt.subplots(2, 1, figsize=(10, 12),
                            gridspec_kw={'height_ratios': [3, 1]},
                            sharex=True)
                            
    ax_pcolormesh = axs[0]  # Upper plot
    ax_profile = axs[1]     # Lower profile plot
    
    fig.subplots_adjust(right=0.8)
    cax = fig.add_axes([0.85, 0.36, 0.05, 0.51])
    
    pcm = ax_pcolormesh.pcolormesh(X, Y, masked_vals, shading="auto", cmap="plasma", norm=norm_method)
    ax_pcolormesh.set_ylabel(y_label, fontsize=18)
    ax_profile.set_xticks(plotting_centers)
    #ax_pcolormesh.set_xticks(np.arange(len(hist_names)) + 0.5)
    
    fig.colorbar(pcm, cax=cax, orientation='vertical')
    
    ax_profile.errorbar(plotting_centers, profile, yerr=asymmetric_errors, fmt='-o', color='#5790fc', capsize=5)
    ax_profile.set_xlabel("HCAL Depth Level")
    ax_profile.set_ylabel("Profile", fontsize=18)
    ax_profile.set_xticklabels(hist_names)
    ax_profile.grid(True, linestyle="--", alpha=0.5)
    
    fig.savefig(outdir + plot_filename + ".png")
    plt.close(fig)
    
###################################################################################################
# MAIN METHOD

def main():

    args = parseArgs()
    
    todaysDate = datetime.date.today().strftime('%Y%m%d')
    
    f = uproot.open(args.input)
    
    hep.style.use("CMS")
    
    pt_min_vals = [0,50,100,150,200,250]
    pt_max_vals = [50,100,150,200,250,300]
    pt_min_vals = []
    pt_max_vals = []

    regions = ["HE", "HB"]
    
    for region in regions:
    
        if region=="HE": range_end = 8
        else: range_end = 5
   
        
        # create leading jet average hcal energy fraction plot
        histograms = [f[f"leadJet_{region}_avgConstituentHcalDepthEF{i}"] for i in range(1,range_end)]
        print(histograms[0])
        plotHCalDepthEFs(histograms,
                         "Lead Jet Average HCAL Energy Fraction",
                         args.outdir,
                         f"leadJet_{region}_avgConstituentHcalDepthEF_{args.dataset}_allPt_{todaysDate}",
                         region)
        for i in range(len(pt_min_vals)):
            plotHCalDepthEFs(histograms,
                             "Lead Jet Average HCAL Energy Fraction",
                             args.outdir,
                             f"leadJet_{region}_avgConstituentHcalDepthEF_{args.dataset}_Pt{pt_min_vals[i]}to{pt_max_vals[i]}_{todaysDate}",
                             region,
                             pt_min=pt_min_vals[i],
                             pt_max=pt_max_vals[i])
            
                         
        # create leading jet pt-weighted average hcal energy fraction plot
        histograms = [f[f"leadJet_{region}_pTWeightedAvgConstituentHcalDepthEF{i}"] for i in range(1,range_end)]
        plotHCalDepthEFs(histograms,
                         "Lead Jet pT-Weighted Average HCAL Energy Fraction",
                         args.outdir,
                         f"leadJet_{region}_pTWeightedAvgConstituentHcalDepthEF_{args.dataset}_allPt_{todaysDate}",
                         region)
        for i in range(len(pt_min_vals)):
            plotHCalDepthEFs(histograms,
                             "Lead Jet pT-Weighted Average HCAL Energy Fraction",
                             args.outdir,
                             f"leadJet_{region}_pTWeightedAvgConstituentHcalDepthEF_{args.dataset}_Pt{pt_min_vals[i]}to{pt_max_vals[i]}_{todaysDate}",
                             region,
                             pt_min=pt_min_vals[i],
                             pt_max=pt_max_vals[i])
                         
        # create leading jet median hcal energy fraction plot
        histograms = [f[f"leadJet_{region}_medConstituentHcalDepthEF{i}"] for i in range(1,range_end)]
        plotHCalDepthEFs(histograms,
                         "Lead Jet Median HCAL Energy Fraction",
                         args.outdir,
                         f"leadJet_{region}_medConstituentHcalDepthEF_{args.dataset}_allPt_{todaysDate}",
                         region)
        for i in range(len(pt_min_vals)):
            plotHCalDepthEFs(histograms,
                             "Lead Jet Median HCAL Energy Fraction",
                             args.outdir,
                             f"leadJet_{region}_medConstituentHcalDepthEF_{args.dataset}_Pt{pt_min_vals[i]}to{pt_max_vals[i]}_{todaysDate}",
                             region,
                             pt_min=pt_min_vals[i],
                             pt_max=pt_max_vals[i])
                         
        # create leading jet minimum hcal energy fraction plot
        histograms = [f[f"leadJet_{region}_minConstituentHcalDepthEF{i}"] for i in range(1,range_end)]
        plotHCalDepthEFs(histograms,
                         "Lead Jet Minimum HCAL Energy Fraction",
                         args.outdir,
                         f"leadJet_{region}_minConstituentHcalDepthEF_{args.dataset}_allPt_{todaysDate}",
                         region)
        for i in range(len(pt_min_vals)):
            plotHCalDepthEFs(histograms,
                             "Lead Jet Minimum HCAL Energy Fraction",
                             args.outdir,
                             f"leadJet_{region}_minConstituentHcalDepthEF_{args.dataset}_Pt{pt_min_vals[i]}to{pt_max_vals[i]}_{todaysDate}",
                             region,
                             pt_min=pt_min_vals[i],
                             pt_max=pt_max_vals[i])
                         
        # create leading jet maximum hcal energy fraction plot
        histograms = [f[f"leadJet_{region}_maxConstituentHcalDepthEF{i}"] for i in range(1,range_end)]
        plotHCalDepthEFs(histograms,
                         "Lead Jet Maximum HCAL Energy Fraction",
                         args.outdir,
                         f"leadJet_{region}_maxConstituentHcalDepthEF_{args.dataset}_allPt_{todaysDate}", 
                         region)
        for i in range(len(pt_min_vals)):
            plotHCalDepthEFs(histograms,
                             "Lead Jet Maximum HCAL Energy Fraction",
                             args.outdir,
                             f"leadJet_{region}_maxConstituentHcalDepthEF_{args.dataset}_Pt{pt_min_vals[i]}to{pt_max_vals[i]}_{todaysDate}",
                             region,
                             pt_min=pt_min_vals[i],
                             pt_max=pt_max_vals[i])
    

###################################################################################################
# RUN SCRIPT

if __name__ == "__main__":
    main()
