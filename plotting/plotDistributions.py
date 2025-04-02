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
    
def plotHcalDepthEFDistributions(sampleDict, histograms, samples, x_label, outdir, plot_filename,
    norm=True, pt_min=None, pt_max=None):

    colors = ['black', 'crimson', 'slateblue', 'tomato', 'seagreen', 'teal', 'darkmagenta']
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    
    max_y = 0.0
    
    for i in range(len(histograms)):
    
        hist_values = histograms[i].values()
        bin_edges = histograms[i].axes[0].edges()
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
        if pt_min is not None and pt_max is not None:
            pt_mask = (histograms[i].axes[1].centers >= pt_min) & (histograms[i].axes[1].centers <= pt_max)
            values = np.sum(hist_values[:, pt_mask], axis=1)
        else:
            values = np.sum(hist_values, axis=1)
            
        errors = np.sqrt(values)
            
        if norm:
            total = np.sum(values)
            values = values / total if total > 0 else values
            errors = errors / total if total > 0 else errors
            
        max_y_current = max(values)
        if max_y < max_y_current: max_y = max_y_current
        
        if colors[i % len(colors)] == 'black':
            linewidth = 2
        else:
            linewidth = 1
            
        ax.errorbar(x=bin_centers, y=values, yerr=errors, linestyle='', color=colors[i % len(colors)])
        ax.errorbar(x=bin_edges[:-1],y=values,drawstyle='steps-post', linewidth=linewidth, label=samples[i],color=colors[i % len(colors)])
        
        
    ax.set_xlabel(x_label)
    ax.set_ylabel("Events [A.U.]" if norm else "Events")
    ax.set_ylim([0,1.6*max_y])
    handles, labels = ax.get_legend_handles_labels()
    unique_labels = {}
    for handle, label in zip(handles, labels):
        unique_labels[label] = handle  # Keep only the last occurrence
    ax.legend(unique_labels.values(), unique_labels.keys(),frameon=False,loc='upper center', fontsize=12,)
    
    fig.savefig(f"{outdir}/{plot_filename}")
    plt.close(fig)

    return
        
    
###################################################################################################
# MAIN METHOD

def main():

    args = parseArgs()
    todaysDate = datetime.date.today().strftime('%Y%m%d')
    hep.style.use("CMS")
    
    sampleDict = constructSampleDict()
    
    pt_min_vals = [0,50,100,150,200,250]
    pt_max_vals = [50,100,150,200,250,300]
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
                    
                    for i in range(1, range_end):
                    
                        histograms = [files[j][f"leadJet_{region}_{category}ConstituentHcalDepthEF{i}"]
                                      for j in range(len(samples))]
                                      
                        plotHcalDepthEFDistributions(
                            sampleDict,
                            histograms,
                            samples,
                            f"Lead Jet {category_names[k]} HCAL Energy Fraction (Depth  {i})",
                            args.outdir,
                            f"leadJet_{region}_{category}ConstituentHCalDepthEF{i}_mDark{mDark_val}_ctau{ctau_val}_{todaysDate}"
                        )


###################################################################################################
# RUN SCRIPT

if __name__ == "__main__":
    main()
