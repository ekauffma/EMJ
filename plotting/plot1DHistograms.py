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
import ROOT

###################################################################################################
# GLOBAL VARIABLES

# uncomment those you want plotted (make sure they have associated colors and names in the subsequent dictionaries)
inputFileDict = {
    "QCD_Bin-Pt-600to800": "histograms/histograms_QCD_Bin-Pt-600to800.root",
#    "EMJ_s-channel_mMed-100_mDark-10_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-10_ctau-1000_unflavored-down.root",
#    "EMJ_s-channel_mMed-100_mDark-10_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-10_ctau-100_unflavored-down.root",
#    "EMJ_s-channel_mMed-100_mDark-10_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-10_ctau-1_unflavored-down.root",
#    "EMJ_s-channel_mMed-100_mDark-20_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-20_ctau-1000_unflavored-down.root",
#    "EMJ_s-channel_mMed-100_mDark-20_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-20_ctau-100_unflavored-down.root",
#    "EMJ_s-channel_mMed-100_mDark-20_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-20_ctau-1_unflavored-down.root",
    "EMJ_s-channel_mMed-250_mDark-10_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-10_ctau-1000_unflavored-down.root",
    "EMJ_s-channel_mMed-250_mDark-10_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-10_ctau-100_unflavored-down.root",
    "EMJ_s-channel_mMed-250_mDark-10_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-10_ctau-1_unflavored-down.root",
    "EMJ_s-channel_mMed-250_mDark-20_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-20_ctau-1000_unflavored-down.root",
    "EMJ_s-channel_mMed-250_mDark-20_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-20_ctau-100_unflavored-down.root",
    "EMJ_s-channel_mMed-250_mDark-20_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-20_ctau-1_unflavored-down.root",
#    "EMJ_s-channel_mMed-500_mDark-10_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-500_mDark-10_ctau-100_unflavored-down.root",
#    "EMJ_s-channel_mMed-500_mDark-10_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-500_mDark-10_ctau-1_unflavored-down.root",
}

datasetNameDict = {
    "QCD_Bin-Pt-600to800": "QCD (pT=600 to 800 GeV)",
    "EMJ_s-channel_mMed-100_mDark-10_ctau-1000_unflavored-down": "EMJ (mMed=100, mDark=10, ctau=1000)",
    "EMJ_s-channel_mMed-100_mDark-10_ctau-100_unflavored-down": "EMJ (mMed=100, mDark=10, ctau=100)",
    "EMJ_s-channel_mMed-100_mDark-10_ctau-1_unflavored-down": "EMJ (mMed=100, mDark=10, ctau=1)",
    "EMJ_s-channel_mMed-100_mDark-20_ctau-1000_unflavored-down": "EMJ (mMed=100, mDark=20, ctau=1000)",
    "EMJ_s-channel_mMed-100_mDark-20_ctau-100_unflavored-down": "EMJ (mMed=100, mDark=20, ctau=100)",
    "EMJ_s-channel_mMed-100_mDark-20_ctau-1_unflavored-down": "EMJ (mMed=100, mDark=20, ctau=1)",
    "EMJ_s-channel_mMed-250_mDark-10_ctau-1000_unflavored-down": "EMJ (mMed=250, mDark=10, ctau=1000)",
    "EMJ_s-channel_mMed-250_mDark-10_ctau-100_unflavored-down": "EMJ (mMed=250, mDark=10, ctau=10)",
    "EMJ_s-channel_mMed-250_mDark-10_ctau-1_unflavored-down": "EMJ (mMed=250, mDark=10, ctau=1)",
    "EMJ_s-channel_mMed-250_mDark-20_ctau-1000_unflavored-down": "EMJ (mMed=250, mDark=20, ctau=1000)",
    "EMJ_s-channel_mMed-250_mDark-20_ctau-100_unflavored-down": "EMJ (mMed=250, mDark=20, ctau=100)",
    "EMJ_s-channel_mMed-250_mDark-20_ctau-1_unflavored-down": "EMJ (mMed=250, mDark=20, ctau=1)",
    "EMJ_s-channel_mMed-500_mDark-10_ctau-100_unflavored-down": "EMJ (mMed=500, mDark=10, ctau=100)", 
    "EMJ_s-channel_mMed-500_mDark-10_ctau-1_unflavored-down": "EMJ (mMed=500, mDark=10, ctau=1)",
}

plotColorDict = {
    "QCD_Bin-Pt-600to800": 1,
    "EMJ_s-channel_mMed-250_mDark-10_ctau-1000_unflavored-down": 50,
    "EMJ_s-channel_mMed-250_mDark-10_ctau-100_unflavored-down": 9,
    "EMJ_s-channel_mMed-250_mDark-10_ctau-1_unflavored-down": 8,
    "EMJ_s-channel_mMed-250_mDark-20_ctau-1000_unflavored-down": 63,
    "EMJ_s-channel_mMed-250_mDark-20_ctau-100_unflavored-down": 94,
    "EMJ_s-channel_mMed-250_mDark-20_ctau-1_unflavored-down": 51,
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
   
def plot1DComparison(histDict, xLabel, xMin, xMax, outFile):

    canvas = ROOT.TCanvas("canvas", "Canvas with TPads", 1000, 600)

    # Define two TPads
    pad1 = ROOT.TPad("pad1", "Pad 1", 0.0, 0.0, 0.8, 1.0)
    pad2 = ROOT.TPad("pad2", "Pad 2", 0.75, 0.0, 1.0, 1.0)

    # Draw the TPads on the canvas
    pad1.Draw()
    pad2.Draw()

    pad1.cd()

    first=True
    histograms = []
    for dataset, hist in histDict.items():
        # normalize
        hist.Scale(1/hist.Integral(0,hist.GetNbinsX()+1))
        
        if first:
            hist.GetXaxis().SetRangeUser(xMin, xMax)
            hist.SetStats(0)
            hist.SetTitle("")
            hist.GetYaxis().SetTitle("Events [A.U.]")
            hist.GetXaxis().SetTitle(xLabel)

        hist.SetMarkerColor(plotColorDict[dataset])
        hist.SetMarkerStyle(20)
        hist.SetLineColor(plotColorDict[dataset])

        histograms.append(hist)

    for i in range(len(histograms)):
        if i==0: histograms[i].Draw("e")
        else: histograms[i].Draw("e same")

    canvas.cd()
    pad2.cd()
    legend = ROOT.TLegend(0, 0.4, 0.9, 0.9)
    dataset_list = list(histDict.keys())
    for i in range(len(histograms)):
        legend.AddEntry(histograms[i], datasetNameDict[dataset_list[i]], "lp")
    legend.SetTextSize(0.055)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.Draw()

    canvas.SaveAs(outFile)

def getHistDict(openFileDict, histName):
    histDict = {}
    for dataset, openFile in openFileDict.items():
        histDict[dataset] = openFile.Get(histName)

    return histDict

    
###################################################################################################
# MAIN METHOD

def main():

    args = parseArgs()
    
    todaysDate = datetime.date.today().strftime('%Y%m%d')

    openFileDict = {}
    for dataset, file in inputFileDict.items():
        openFileDict[dataset] = ROOT.TFile(file)
   
    # lead jet pt plot
    histDict = getHistDict(openFileDict, "leadJet_Pt")
    plot1DComparison(histDict, r"Lead Jet p_T [GeV]", 0, 1000, f"{args.outdir}/leadjetpt_{todaysDate}.pdf")
    
    # lead jet eta plot
    histDict = getHistDict(openFileDict, "leadJet_Eta")
    plot1DComparison(histDict, r"Lead Jet #eta", -3, 3, f"{args.outdir}/leadjeteta_{todaysDate}.pdf")
    
    # lead jet phi plot
    histDict = getHistDict(openFileDict, "leadJet_Phi")
    plot1DComparison(histDict, r"Lead Jet #phi", -3, 3, f"{args.outdir}/leadjetphi_{todaysDate}.pdf")
    
    # lead jet E plot
    histDict = getHistDict(openFileDict, "leadJet_E")
    plot1DComparison(histDict, r"Lead Jet E [GeV]", 0, 1000, f"{args.outdir}/leadjete_{todaysDate}.pdf")
    
    # lead jet nConstituent plot
    histDict = getHistDict(openFileDict, "leadJet_nConstituent")
    plot1DComparison(histDict, r"Lead Jet Number of Constituents", 0, 200, f"{args.outdir}/leadjetnconstituent_{todaysDate}.pdf")
    
    for file in openFileDict.values():
        file.Close()

###################################################################################################
# RUN SCRIPT

if __name__ == "__main__":
    main()
