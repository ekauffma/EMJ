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
    X, Y = np.meshgrid(depth_edges, bin_edges)
    values_2d = np.array([hist.values() for hist in histograms]).T
    
    masked_vals = np.ma.masked_where(values_2d == 0, values_2d)
    
    fig, ax = plt.subplots(figsize=(10, 10))
        
    pcm = ax.pcolormesh(X, Y, masked_vals, shading="auto", cmap="plasma", norm=norm_method)
    
    ax.set_xlabel("HCAL Depth Level")
    ax.set_ylabel(y_label)
    ax.set_xticks(np.arange(len(hist_names)) + 0.5)
    ax.set_xticklabels(hist_names, rotation=45)
    cbar = fig.colorbar(pcm, ax=ax, orientation='vertical')

    fig.savefig(outdir + plot_filename)
    
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
                     f"leadJet_avgConstituentHcalDepthEF_{args.dataset}_{todaysDate}.pdf")
                     
    # create leading jet pt-weighted average hcal energy fraction plot
    histograms = [f[f"leadJet_pTWeightedAvgConstituentHcalDepthEF{i}"] for i in range(1,8)]
    plotHCalDepthEFs(histograms,
                     "Lead Jet pT-Weighted Average HCAL Energy Fraction",
                     args.outdir,
                     f"leadJet_pTWeightedAvgConstituentHcalDepthEF_{args.dataset}_{todaysDate}.pdf")
                     
    # create leading jet median hcal energy fraction plot
    histograms = [f[f"leadJet_medConstituentHcalDepthEF{i}"] for i in range(1,8)]
    plotHCalDepthEFs(histograms,
                     "Lead Jet Median HCAL Energy Fraction",
                     args.outdir,
                     f"leadJet_medConstituentHcalDepthEF_{args.dataset}_{todaysDate}.pdf")
                     
    # create leading jet minimum hcal energy fraction plot
    histograms = [f[f"leadJet_minConstituentHcalDepthEF{i}"] for i in range(1,8)]
    plotHCalDepthEFs(histograms,
                     "Lead Jet Minimum HCAL Energy Fraction",
                     args.outdir,
                     f"leadJet_minConstituentHcalDepthEF_{args.dataset}_{todaysDate}.pdf")
                     
    # create leading jet maximum hcal energy fraction plot
    histograms = [f[f"leadJet_maxConstituentHcalDepthEF{i}"] for i in range(1,8)]
    plotHCalDepthEFs(histograms,
                     "Lead Jet Maximum HCAL Energy Fraction",
                     args.outdir,
                     f"leadJet_maxConstituentHcalDepthEF_{args.dataset}_{todaysDate}.pdf")
    

###################################################################################################
# RUN SCRIPT

if __name__ == "__main__":
    main()
