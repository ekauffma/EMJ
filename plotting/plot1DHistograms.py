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
    "QCD_Combined": "histograms/histograms_QCD_Combined.root",
#    "QCD_Bin-Pt-30to50": "histograms/histograms_QCD_Bin-Pt-30to50.root",
#    "QCD_Bin-Pt-50to80": "histograms/histograms_QCD_Bin-Pt-50to80.root",
#    "QCD_Bin-Pt-80to120": "histograms/histograms_QCD_Bin-Pt-80to120.root",
#    "QCD_Bin-Pt-120to170": "histograms/histograms_QCD_Bin-Pt-120to170.root",
#    "QCD_Bin-Pt-170to300": "histograms/histograms_QCD_Bin-Pt-170to300.root",
#    "QCD_Bin-Pt-300to470": "histograms/histograms_QCD_Bin-Pt-300to470.root",
#    "QCD_Bin-Pt-470to600": "histograms/histograms_QCD_Bin-Pt-470to600.root",
#    "QCD_Bin-Pt-600to800": "histograms/histograms_QCD_Bin-Pt-600to800.root",
#    "QCD_Bin-Pt-800to1000": "histograms/histograms_QCD_Bin-Pt-800to1000.root",
#    "QCD_Bin-Pt-1000": "histograms/histograms_QCD_Bin-Pt-1000.root",
#    "EMJ_s-channel_mMed-100_mDark-10_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-10_ctau-1000_unflavored-down.root",
#    "EMJ_s-channel_mMed-100_mDark-10_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-10_ctau-100_unflavored-down.root",
#    "EMJ_s-channel_mMed-100_mDark-10_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-10_ctau-1_unflavored-down.root",
#    "EMJ_s-channel_mMed-100_mDark-20_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-20_ctau-1000_unflavored-down.root",
#    "EMJ_s-channel_mMed-100_mDark-20_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-20_ctau-100_unflavored-down.root",
#    "EMJ_s-channel_mMed-100_mDark-20_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-100_mDark-20_ctau-1_unflavored-down.root",
#    "EMJ_s-channel_mMed-250_mDark-10_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-10_ctau-1000_unflavored-down.root",
#    "EMJ_s-channel_mMed-250_mDark-10_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-10_ctau-100_unflavored-down.root",
#    "EMJ_s-channel_mMed-250_mDark-10_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-10_ctau-1_unflavored-down.root",
#    "EMJ_s-channel_mMed-250_mDark-20_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-20_ctau-1000_unflavored-down.root",
#    "EMJ_s-channel_mMed-250_mDark-20_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-20_ctau-100_unflavored-down.root",
#    "EMJ_s-channel_mMed-250_mDark-20_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-250_mDark-20_ctau-1_unflavored-down.root",
#    "EMJ_s-channel_mMed-500_mDark-10_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-500_mDark-10_ctau-1000_unflavored-down.root",
#    "EMJ_s-channel_mMed-500_mDark-10_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-500_mDark-10_ctau-100_unflavored-down.root",
#    "EMJ_s-channel_mMed-500_mDark-10_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-500_mDark-10_ctau-1_unflavored-down.root",
#    "EMJ_s-channel_mMed-500_mDark-20_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-500_mDark-20_ctau-1000_unflavored-down.root",
#    "EMJ_s-channel_mMed-500_mDark-20_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-500_mDark-20_ctau-100_unflavored-down.root",
#    "EMJ_s-channel_mMed-500_mDark-20_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-500_mDark-20_ctau-1_unflavored-down.root",
#    "EMJ_s-channel_mMed-750_mDark-10_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-750_mDark-10_ctau-1000_unflavored-down.root",
#    "EMJ_s-channel_mMed-750_mDark-10_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-750_mDark-10_ctau-100_unflavored-down.root",
#    "EMJ_s-channel_mMed-750_mDark-10_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-750_mDark-10_ctau-1_unflavored-down.root",
#    "EMJ_s-channel_mMed-750_mDark-20_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-750_mDark-20_ctau-1000_unflavored-down.root",
#    "EMJ_s-channel_mMed-750_mDark-20_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-750_mDark-20_ctau-100_unflavored-down.root",
#    "EMJ_s-channel_mMed-750_mDark-20_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-750_mDark-20_ctau-1_unflavored-down.root",
#    "EMJ_s-channel_mMed-1000_mDark-10_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-1000_mDark-10_ctau-1000_unflavored-down.root",
#    "EMJ_s-channel_mMed-1000_mDark-10_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-1000_mDark-10_ctau-100_unflavored-down.root",
#    "EMJ_s-channel_mMed-1000_mDark-10_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-1000_mDark-10_ctau-1_unflavored-down.root",
#    "EMJ_s-channel_mMed-1000_mDark-20_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-1000_mDark-20_ctau-1000_unflavored-down.root",
#    "EMJ_s-channel_mMed-1000_mDark-20_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-1000_mDark-20_ctau-100_unflavored-down.root",
#    "EMJ_s-channel_mMed-1000_mDark-20_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-1000_mDark-20_ctau-1_unflavored-down.root",
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-1500_mDark-10_ctau-1000_unflavored-down.root",
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-1500_mDark-10_ctau-100_unflavored-down.root",
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-1500_mDark-10_ctau-1_unflavored-down.root",
    "EMJ_s-channel_mMed-1500_mDark-20_ctau-1000_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-1500_mDark-20_ctau-1000_unflavored-down.root",
    "EMJ_s-channel_mMed-1500_mDark-20_ctau-100_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-1500_mDark-20_ctau-100_unflavored-down.root",
    "EMJ_s-channel_mMed-1500_mDark-20_ctau-1_unflavored-down": "histograms/histograms_EMJ_s-channel_mMed-1500_mDark-20_ctau-1_unflavored-down.root",
}

datasetNameDict = {
    "QCD_Combined": "QCD",
    "QCD_Bin-Pt-30to50": "QCD (pT=30 to 50 GeV)",
    "QCD_Bin-Pt-50to80": "QCD (pT=50 to 80 GeV)",
    "QCD_Bin-Pt-80to120": "QCD (pT=80 to 120 GeV)",
    "QCD_Bin-Pt-120to170": "QCD (pT=120 to 170 GeV)",
    "QCD_Bin-Pt-170to300": "QCD (pT=170 to 300 GeV)",
    "QCD_Bin-Pt-300to470": "QCD (pT=300 to 470 GeV)",
    "QCD_Bin-Pt-470to600": "QCD (pT=470 to 600 GeV)",
    "QCD_Bin-Pt-600to800": "QCD (pT=600 to 800 GeV)",
    "QCD_Bin-Pt-800to1000": "QCD (pT=800 to 1000 GeV)",
    "QCD_Bin-Pt-1000": "QCD (pT>1000 GeV)",
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
    "EMJ_s-channel_mMed-500_mDark-10_ctau-1000_unflavored-down": "EMJ (mMed=500, mDark=10, ctau=1000)",
    "EMJ_s-channel_mMed-500_mDark-10_ctau-100_unflavored-down": "EMJ (mMed=500, mDark=10, ctau=10)",
    "EMJ_s-channel_mMed-500_mDark-10_ctau-1_unflavored-down": "EMJ (mMed=500, mDark=10, ctau=1)",
    "EMJ_s-channel_mMed-500_mDark-20_ctau-1000_unflavored-down": "EMJ (mMed=500, mDark=20, ctau=1000)",
    "EMJ_s-channel_mMed-500_mDark-20_ctau-100_unflavored-down": "EMJ (mMed=500, mDark=20, ctau=100)",
    "EMJ_s-channel_mMed-500_mDark-20_ctau-1_unflavored-down": "EMJ (mMed=250, mDark=20, ctau=1)",
    "EMJ_s-channel_mMed-750_mDark-10_ctau-1000_unflavored-down": "EMJ (mMed=750, mDark=10, ctau=1000)",
    "EMJ_s-channel_mMed-750_mDark-10_ctau-100_unflavored-down": "EMJ (mMed=750, mDark=10, ctau=10)",
    "EMJ_s-channel_mMed-750_mDark-10_ctau-1_unflavored-down": "EMJ (mMed=750, mDark=10, ctau=1)",
    "EMJ_s-channel_mMed-750_mDark-20_ctau-1000_unflavored-down": "EMJ (mMed=750, mDark=20, ctau=1000)",
    "EMJ_s-channel_mMed-750_mDark-20_ctau-100_unflavored-down": "EMJ (mMed=750, mDark=20, ctau=100)",
    "EMJ_s-channel_mMed-750_mDark-20_ctau-1_unflavored-down": "EMJ (mMed=750, mDark=20, ctau=1)",
    "EMJ_s-channel_mMed-1000_mDark-10_ctau-1000_unflavored-down": "EMJ (mMed=1000, mDark=10, ctau=1000)",
    "EMJ_s-channel_mMed-1000_mDark-10_ctau-100_unflavored-down": "EMJ (mMed=1000, mDark=10, ctau=10)",
    "EMJ_s-channel_mMed-1000_mDark-10_ctau-1_unflavored-down": "EMJ (mMed=1000, mDark=10, ctau=1)",
    "EMJ_s-channel_mMed-1000_mDark-20_ctau-1000_unflavored-down": "EMJ (mMed=1000, mDark=20, ctau=1000)",
    "EMJ_s-channel_mMed-1000_mDark-20_ctau-100_unflavored-down": "EMJ (mMed=1000, mDark=20, ctau=100)",
    "EMJ_s-channel_mMed-1000_mDark-20_ctau-1_unflavored-down": "EMJ (mMed=1000, mDark=20, ctau=1)",
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-1000_unflavored-down": "EMJ (mMed=1500, mDark=10, ctau=1000)",
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-100_unflavored-down": "EMJ (mMed=1500, mDark=10, ctau=10)",
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-1_unflavored-down": "EMJ (mMed=1500, mDark=10, ctau=1)",
    "EMJ_s-channel_mMed-1500_mDark-20_ctau-1000_unflavored-down": "EMJ (mMed=1500, mDark=20, ctau=1000)",
    "EMJ_s-channel_mMed-1500_mDark-20_ctau-100_unflavored-down": "EMJ (mMed=1500, mDark=20, ctau=100)",
    "EMJ_s-channel_mMed-1500_mDark-20_ctau-1_unflavored-down": "EMJ (mMed=1500, mDark=20, ctau=1)",
}

plotColorDict = {
    "QCD_Combined": 1,
    "QCD_Bin-Pt-30to50": 50,
    "QCD_Bin-Pt-50to80": 9,
    "QCD_Bin-Pt-80to120": 8,
    "QCD_Bin-Pt-120to170": 63,
    "QCD_Bin-Pt-170to300": 94,
    "QCD_Bin-Pt-300to470": 51,
    "QCD_Bin-Pt-470to600": 38,
    "QCD_Bin-Pt-600to800": 67,
    "QCD_Bin-Pt-800to1000": 6,
    "QCD_Bin-Pt-1000": 71,
    "EMJ_s-channel_mMed-100_mDark-10_ctau-1000_unflavored-down": 50,
    "EMJ_s-channel_mMed-100_mDark-10_ctau-100_unflavored-down": 9,
    "EMJ_s-channel_mMed-100_mDark-10_ctau-1_unflavored-down": 8,
    "EMJ_s-channel_mMed-100_mDark-20_ctau-1000_unflavored-down": 63,
    "EMJ_s-channel_mMed-100_mDark-20_ctau-100_unflavored-down": 94,
    "EMJ_s-channel_mMed-100_mDark-20_ctau-1_unflavored-down": 51,
    "EMJ_s-channel_mMed-250_mDark-10_ctau-1000_unflavored-down": 50,
    "EMJ_s-channel_mMed-250_mDark-10_ctau-100_unflavored-down": 9,
    "EMJ_s-channel_mMed-250_mDark-10_ctau-1_unflavored-down": 8,
    "EMJ_s-channel_mMed-250_mDark-20_ctau-1000_unflavored-down": 63,
    "EMJ_s-channel_mMed-250_mDark-20_ctau-100_unflavored-down": 94,
    "EMJ_s-channel_mMed-250_mDark-20_ctau-1_unflavored-down": 51,
    "EMJ_s-channel_mMed-500_mDark-10_ctau-1000_unflavored-down": 50,
    "EMJ_s-channel_mMed-500_mDark-10_ctau-100_unflavored-down": 9,
    "EMJ_s-channel_mMed-500_mDark-10_ctau-1_unflavored-down": 8,
    "EMJ_s-channel_mMed-500_mDark-20_ctau-1000_unflavored-down": 63,
    "EMJ_s-channel_mMed-500_mDark-20_ctau-100_unflavored-down": 94,
    "EMJ_s-channel_mMed-500_mDark-20_ctau-1_unflavored-down": 51,
    "EMJ_s-channel_mMed-750_mDark-10_ctau-1000_unflavored-down": 50,
    "EMJ_s-channel_mMed-750_mDark-10_ctau-100_unflavored-down": 9,
    "EMJ_s-channel_mMed-750_mDark-10_ctau-1_unflavored-down": 8,
    "EMJ_s-channel_mMed-750_mDark-20_ctau-1000_unflavored-down": 63,
    "EMJ_s-channel_mMed-750_mDark-20_ctau-100_unflavored-down": 94,
    "EMJ_s-channel_mMed-750_mDark-20_ctau-1_unflavored-down": 51,
    "EMJ_s-channel_mMed-1000_mDark-10_ctau-1000_unflavored-down": 50,
    "EMJ_s-channel_mMed-1000_mDark-10_ctau-100_unflavored-down": 9,
    "EMJ_s-channel_mMed-1000_mDark-10_ctau-1_unflavored-down": 8,
    "EMJ_s-channel_mMed-1000_mDark-20_ctau-1000_unflavored-down": 63,
    "EMJ_s-channel_mMed-1000_mDark-20_ctau-100_unflavored-down": 94,
    "EMJ_s-channel_mMed-1000_mDark-20_ctau-1_unflavored-down": 51,
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-1000_unflavored-down": 50,
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-100_unflavored-down": 9,
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-1_unflavored-down": 8,
    "EMJ_s-channel_mMed-1500_mDark-20_ctau-1000_unflavored-down": 63,
    "EMJ_s-channel_mMed-1500_mDark-20_ctau-100_unflavored-down": 94,
    "EMJ_s-channel_mMed-1500_mDark-20_ctau-1_unflavored-down": 51,
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
    y_max = 0.0
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
        hist.SetMarkerStyle(6)
        hist.SetLineColor(plotColorDict[dataset])
        
        y_max_current = hist.GetBinContent(hist.GetMaximumBin())
        if y_max_current>y_max:
            y_max = y_max_current

        histograms.append(hist)

    histograms[0].GetYaxis().SetRangeUser(0, 1.15*y_max)
    for i in range(len(histograms)):
        if i==0: histograms[i].Draw("hist")
        else: histograms[i].Draw("hist same")
        histograms[i].Draw("e same")

    canvas.cd()
    pad2.cd()
    legend = ROOT.TLegend(0, 0.35, 0.9, 0.9)
    dataset_list = list(histDict.keys())
    for i in range(len(histograms)):
        legend.AddEntry(histograms[i], datasetNameDict[dataset_list[i]], "lp")
    legend.SetTextSize(0.055)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.Draw()

    # canvas.SaveAs(outFile + ".pdf")
    canvas.SaveAs(outFile + ".png")

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
    plot1DComparison(histDict, r"Lead Jet $p_T [GeV]$", 0, 800, f"{args.outdir}/leadjetpt_{todaysDate}")
    
    # lead jet eta plot
    histDict = getHistDict(openFileDict, "leadJet_Eta")
    plot1DComparison(histDict, r"Lead Jet $\eta$", -3, 3, f"{args.outdir}/leadjeteta_{todaysDate}")
    
    # lead jet phi plot
    histDict = getHistDict(openFileDict, "leadJet_Phi")
    plot1DComparison(histDict, r"Lead Jet $\phi$", -3, 3, f"{args.outdir}/leadjetphi_{todaysDate}")
    
    # lead jet E plot
    histDict = getHistDict(openFileDict, "leadJet_E")
    plot1DComparison(histDict, r"Lead Jet $E$ [GeV]", 0, 1000, f"{args.outdir}/leadjete_{todaysDate}")
    
    # lead jet nConstituent plot
    histDict = getHistDict(openFileDict, "leadJet_nConstituent")
    plot1DComparison(histDict, r"Lead Jet Number of Constituents", 0, 100, f"{args.outdir}/leadjetnconstituent_{todaysDate}")
    
    for file in openFileDict.values():
        file.Close()

###################################################################################################
# RUN SCRIPT

if __name__ == "__main__":
    main()
