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
   
def plot1DComparison(sampleDict, histName, xLabel, xMin, xMax, outFile, log=False):

    openFileDict = {}
    for dataset, v in sampleDict.items():
        file = v["histFile"]
        openFileDict[dataset] = ROOT.TFile(file)

    histDict = getHistDict(openFileDict, histName)
    colors = [64, 95, 50, 51, 6, 8, 9]

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
    i = 0
    for dataset, hist in histDict.items():
        # normalize
        hist.Scale(1/hist.Integral(0,hist.GetNbinsX()+1))
        
        if first:
            hist.GetXaxis().SetRangeUser(xMin, xMax)
            hist.SetStats(0)
            hist.SetTitle("")
            hist.GetYaxis().SetTitle("Events [A.U.]")
            hist.GetXaxis().SetTitle(xLabel)
            hist.SetLineWidth(2)
            first = False

        print(dataset)
        if dataset=="QCD_Combined":
            hist.SetMarkerColor(1)
            hist.SetLineColor(1)
            hist.SetLineWidth(2)
        else:
            hist.SetMarkerColor(colors[i])
            hist.SetLineColor(colors[i])
       
        hist.SetMarkerStyle(6)

        y_max_current = hist.GetBinContent(hist.GetMaximumBin())
        if y_max_current>y_max:
            y_max = y_max_current

        histograms.append(hist)
        i = i+1

    if log:
        histograms[0].GetYaxis().SetRangeUser(1e-3, 10*y_max)
    else:
        histograms[0].GetYaxis().SetRangeUser(0, 1.15*y_max)

    for i in range(len(histograms)):
        if i==0: histograms[i].Draw("hist")
        else: histograms[i].Draw("hist same")
        histograms[i].Draw("e same")

    if log:
        pad1.SetLogy()

    canvas.cd()
    pad2.cd()
    legend = ROOT.TLegend(0, 0.35, 0.9, 0.9)
    dataset_list = list(histDict.keys())
    for i in range(len(histograms)):
        legend.AddEntry(histograms[i], sampleDict[dataset_list[i]]["printName"], "lp")
    legend.SetTextSize(0.055)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.Draw()

    # canvas.SaveAs(outFile + ".pdf")
    canvas.SaveAs(outFile + ".png")

    del histDict
    del canvas

    for file in openFileDict.values():
        file.Close()

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

    sampleDict = constructSampleDict()

    ctau_list = [2000]
    mDark_list = [10, 20]

    for ctau in ctau_list:
        for mDark in mDark_list:
            tempDict = {}
            for k, v in sampleDict.items():
                if f"ctau-{ctau}_" in k:
                    if f"mDark-{mDark}" in k:
                        tempDict[k] = v
            tempDict["QCD_Combined"] = sampleDict["QCD_Combined"]
    
            print(tempDict)
            print()
    
            # lead jet pt plot
            plot1DComparison(tempDict, "leadJet_Pt", r"Lead Jet $p_T [GeV]$", 0, 800, f"{args.outdir}/leadjetpt_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True)
            
            # lead jet eta plot
            plot1DComparison(tempDict, "leadJet_Eta", r"Lead Jet $\eta$", -3, 3, f"{args.outdir}/leadjeteta_mDark-{mDark}_ctau-{ctau}_{todaysDate}")
            
            # lead jet phi plot
            plot1DComparison(tempDict, "leadJet_Phi", r"Lead Jet $\phi$", -3, 3, f"{args.outdir}/leadjetphi_mDark-{mDark}_ctau-{ctau}_{todaysDate}")
            
            # lead jet E plot
            plot1DComparison(tempDict, "leadJet_E", r"Lead Jet $E$ [GeV]", 0, 1000, f"{args.outdir}/leadjete_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True)
            
            # lead jet nConstituent plot
            plot1DComparison(tempDict, "leadJet_nConstituent", r"Lead Jet Number of Constituents", 0, 100, f"{args.outdir}/leadjetnconstituent_mDark-{mDark}_ctau-{ctau}_{todaysDate}")
    

###################################################################################################
# RUN SCRIPT

if __name__ == "__main__":
    main()
