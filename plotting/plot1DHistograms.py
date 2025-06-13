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
   
def plot1DComparison(sampleDict, histName, xLabel, xMin, xMax, outFile, log=False, normStart=0):

    print("Hist Name = ", histName)

    openFileDict = {}
    for dataset, v in sampleDict.items():
        file = v["histFile"]
        openFileDict[dataset] = ROOT.TFile(file)

    histDict = {}
    for k, v in getHistDict(openFileDict, histName).items():
        h = v.Clone(f"{v.GetName()}_{k}")
        h.SetDirectory(0)
        h.Scale(1/h.Integral(normStart,h.GetNbinsX()+1))
        histDict[k] = h
    colors = [64, 95, 50, 51, 6, 8, 9]

    canvas = ROOT.TCanvas("c", "c", 1000, 600)

    pad1 = ROOT.TPad("pad1", "Pad 1", 0.0, 0.0, 0.65, 1.0)
    pad2 = ROOT.TPad("pad2", "Pad 2", 0.65, 0.0, 1.0, 1.0)
    pad1.SetLeftMargin(0.15)
    pad1.SetRightMargin(0.02)  # small right margin
    pad2.SetLeftMargin(0.02)
    pad2.SetRightMargin(0.05)

    # Draw the TPads on the canvas
    pad1.Draw()
    pad2.Draw()

    pad1.cd()
    
    first=True
    histograms = []
    y_max = 0.0
    i = 0
    for dataset, hist in histDict.items():
        if first:
            hist.GetXaxis().SetRangeUser(xMin, xMax)
            hist.SetStats(0)
            hist.SetTitle("")
            hist.GetYaxis().SetTitle("Events [A.U.]")
            hist.GetXaxis().SetTitle(xLabel)
            first = False

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
        histograms[0].GetYaxis().SetRangeUser(1e-5, 2*y_max)
    else:
        histograms[0].GetYaxis().SetRangeUser(0, 1.15*y_max)

    for i in range(len(histograms)):
        if i==0: histograms[i].Draw("hist e")
        else: histograms[i].Draw("hist e same")

    if log:
        pad1.SetLogy()

    canvas.cd()
    pad2.cd()
    legend = ROOT.TLegend(0, 0.6, 0.9, 0.9)
    dataset_list = list(histDict.keys())
    for i in range(len(histograms)):
        legend.AddEntry(histograms[i], sampleDict[dataset_list[i]]["printName"], "lp")
    legend.SetTextSize(0.04)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.Draw()

    # canvas.SaveAs(outFile + ".pdf")
    canvas.Update()
    canvas.SaveAs(outFile + ".png")

    #for hist in histograms:
    #    hist.Delete()
    #del histDict
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

    ctau_list = [1, 100, 1000, 1500, 2000]
    mDark_list = [10]
    ctau_list = [1, 1500]

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
            for region in ["HE", "HB"]:
    
                plot1DComparison(tempDict, f"leadJet_{region}_Pt", "Lead Jet p_{T} [GeV]" + f" ({region})", 0, 800,
                                 f"{args.outdir}/leadjetpt_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True)
                plot1DComparison(tempDict, f"leadJet_{region}_Eta", "Lead Jet #eta" + f" ({region})", -3, 3,
                                 f"{args.outdir}/leadjeteta_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}")
                plot1DComparison(tempDict, f"leadJet_{region}_Phi", "Lead Jet #phi" + f" ({region})", -3, 3,
                                 f"{args.outdir}/leadjetphi_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}")
                plot1DComparison(tempDict, f"leadJet_{region}_E", r"Lead Jet $E$ [GeV]" + f" ({region})", 0, 1000,
                                 f"{args.outdir}/leadjete_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True)
                plot1DComparison(tempDict, f"leadJet_{region}_nConstituent", r"Lead Jet Number of Constituents" + f" ({region})", 0, 100,
                                 f"{args.outdir}/leadjetnconstituent_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}")
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_frac0", "Lead Jet Fraction of RecHits with auxTDC=0" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCfrac0_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_frac1", "Lead Jet Fraction of RecHits with auxTDC=1" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCfrac1_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_frac2", "Lead Jet Fraction of RecHits with auxTDC=2" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCfrac2_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_frac3", "Lead Jet Fraction of RecHits with auxTDC=3" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCfrac3_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_energyWeightedFrac0", "Lead Jet Energy-Weighted Fraction of RecHits with auxTDC=0" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCEnergyWeightedFrac0_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_energyWeightedFrac1", "Lead Jet Energy-Weighted Fraction of RecHits with auxTDC=1" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCEnergyWeightedFrac1_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_energyWeightedFrac2", "Lead Jet Energy-Weighted Fraction of RecHits with auxTDC=2" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCEnergyWeightedFrac2_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_energyWeightedFrac3", "Lead Jet Energy-Weighted Fraction of RecHits with auxTDC=3" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCEnergyWeightedFrac3_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_energy_avg", "Lead Jet Average RecHit Energy [GeV]" + f" ({region})", 0, 1, 
                                 f"{args.outdir}/leadjet_energyavg_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_energy_max", "Lead Jet Maximum RecHit Energy [GeV]" + f" ({region})", 0, 10, 
                                 f"{args.outdir}/leadjet_energymax_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_energy_min", "Lead Jet Minimum RecHit Energy [GeV]" + f" ({region})", 0, 0.5, 
                                 f"{args.outdir}/leadjet_energymin_{region}_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True)
                plot1DComparison(tempDict, f"recHit_{region}_E_max", r"Maximum per-Event RecHit Energy [GeV]" + f" ({region})", 0, 10,
                                 f"{args.outdir}/recHit_{region}_E_max_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True)
                plot1DComparison(tempDict, f"recHit_{region}_E_min", r"Minimum per-Event RecHit Energy [GeV]" + f" ({region})", 0, 0.5,
                                 f"{args.outdir}/recHit_{region}_E_min_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True)
                plot1DComparison(tempDict, f"recHit_{region}_E_median", r"Median per-Event RecHit Energy [GeV]" + f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_E_median_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True)
                plot1DComparison(tempDict, f"recHit_{region}_E_mean", r"Mean per-Event RecHit Energy [GeV]" + f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_E_mean_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_frac0", r"Fraction of RecHit auxTDC 3 equal to 0"+ f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_auxTDC3_frac0_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_frac1", r"Fraction of RecHit auxTDC 3 equal to 1"+ f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_auxTDC3_frac1_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_frac2", r"Fraction of RecHit auxTDC 3 equal to 2"+ f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_auxTDC3_frac2_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_frac3", r"Fraction of RecHit auxTDC 3 equal to 3"+ f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_auxTDC3_frac3_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_energyWeightedFrac0", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 0" + f" ({region})", 0, 1, 
                                 f"{args.outdir}/recHit_{region}_auxTDC3_energyWeightedFrac0_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_energyWeightedFrac1", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 1" + f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_auxTDC3_energyWeightedFrac1_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_energyWeightedFrac2", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 2" + f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_auxTDC3_energyWeightedFrac2_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_energyWeightedFrac3", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 3" + f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_auxTDC3_energyWeightedFrac3_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)
                
            plot1DComparison(tempDict, "leadJet_Pt", "Lead Jet p_{T} [GeV]", 0, 800, f"{args.outdir}/leadjetpt_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True)
            plot1DComparison(tempDict, "leadJet_Eta", "Lead Jet #eta", -3, 3, f"{args.outdir}/leadjeteta_mDark-{mDark}_ctau-{ctau}_{todaysDate}")
            plot1DComparison(tempDict, "leadJet_Phi", "Lead Jet #phi", -3, 3, f"{args.outdir}/leadjetphi_mDark-{mDark}_ctau-{ctau}_{todaysDate}")
            plot1DComparison(tempDict, "leadJet_E", "Lead Jet E [GeV]", 0, 1000, f"{args.outdir}/leadjete_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True)
            plot1DComparison(tempDict, "leadJet_nConstituent", "Lead Jet Number of Constituents", 0, 100, f"{args.outdir}/leadjetnconstituent_mDark-{mDark}_ctau-{ctau}_{todaysDate}")
            plot1DComparison(tempDict, f"leadJet_JetRecHit_auxTDC3_frac0", "Lead Jet Fraction of RecHits with auxTDC=0", 0, 1, f"{args.outdir}/leadjet_auxTDCfrac0_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True, normStart=1)
            plot1DComparison(tempDict, f"leadJet_JetRecHit_auxTDC3_frac1", "Lead Jet Fraction of RecHits with auxTDC=1", 0, 1, f"{args.outdir}/leadjet_auxTDCfrac1_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True, normStart=1)
            plot1DComparison(tempDict, f"leadJet_JetRecHit_auxTDC3_frac2", "Lead Jet Fraction of RecHits with auxTDC=2", 0, 1, f"{args.outdir}/leadjet_auxTDCfrac2_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True, normStart=1)
            plot1DComparison(tempDict, f"leadJet_JetRecHit_auxTDC3_frac3", "Lead Jet Fraction of RecHits with auxTDC=3", 0, 1, f"{args.outdir}/leadjet_auxTDCfrac3_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True, normStart=1)
            plot1DComparison(tempDict, f"leadJet_JetRecHit_energy_avg", "Lead Jet Average RecHit Energy [GeV]", 0, 1, f"{args.outdir}/leadjet_energyavg_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True)
            plot1DComparison(tempDict, f"leadJet_JetRecHit_energy_max", "Lead Jet Maximum RecHit Energy [GeV]", 0, 10, f"{args.outdir}/leadjet_energymax_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True)
            plot1DComparison(tempDict, f"leadJet_JetRecHit_energy_min", "Lead Jet Minimum RecHit Energy [GeV]", 0, 0.5, f"{args.outdir}/leadjet_energymin_mDark-{mDark}_ctau-{ctau}_{todaysDate}", log=True)
            plot1DComparison(tempDict, f"recHit_E_max", r"Maximum per-Event RecHit Energy [GeV]", 0, 10, f"{args.outdir}/recHit_E_max_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True)
            plot1DComparison(tempDict, f"recHit_E_min", r"Minimum per-Event RecHit Energy [GeV]", 0, 0.5, f"{args.outdir}/recHit_E_min_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True)
            plot1DComparison(tempDict, f"recHit_E_median", r"Median per-Event RecHit Energy [GeV]", 0, 1,f"{args.outdir}/recHit_E_median_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True)
            plot1DComparison(tempDict, f"recHit_E_mean", r"Mean per-Event RecHit Energy [GeV]", 0, 1, f"{args.outdir}/recHit_E_mean_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True)
            plot1DComparison(tempDict, f"recHit_auxTDC3_frac0", r"Fraction of RecHit auxTDC 3 equal to 0", 0, 1, f"{args.outdir}/recHit_auxTDC3_frac0_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)
            plot1DComparison(tempDict, f"recHit_auxTDC3_frac1", r"Fraction of RecHit auxTDC 3 equal to 1", 0, 1, f"{args.outdir}/recHit_auxTDC3_frac1_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)
            plot1DComparison(tempDict, f"recHit_auxTDC3_frac2", r"Fraction of RecHit auxTDC 3 equal to 2", 0, 1, f"{args.outdir}/recHit_auxTDC3_frac2_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)
            plot1DComparison(tempDict, f"recHit_auxTDC3_frac3", r"Fraction of RecHit auxTDC 3 equal to 3", 0, 1, f"{args.outdir}/recHit_auxTDC3_frac3_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)
            plot1DComparison(tempDict, f"recHit_auxTDC3_energyWeightedFrac0", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 0", 0, 1, f"{args.outdir}/recHit_auxTDC3_energyWeightedFrac0_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)
            plot1DComparison(tempDict, f"recHit_auxTDC3_energyWeightedFrac1", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 1", 0, 1, f"{args.outdir}/recHit_auxTDC3_energyWeightedFrac1_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)
            plot1DComparison(tempDict, f"recHit_auxTDC3_energyWeightedFrac2", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 2", 0, 1, f"{args.outdir}/recHit_auxTDC3_energyWeightedFrac2_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)
            plot1DComparison(tempDict, f"recHit_auxTDC3_energyWeightedFrac3", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 3", 0, 1, f"{args.outdir}/recHit_auxTDC3_energyWeightedFrac3_mDark-{mDark}_ctau-{ctau}_{todaysDate}",log=True, normStart=1)

    
    mMed_list = [100, 250, 500, 750, 1000, 1500]
    mMed_list = [1500]

    for mMed in mMed_list:
        for mDark in mDark_list:
            tempDict = {}
            for k, v in sampleDict.items():
                if f"mMed-{mMed}_" in k:
                    if f"mDark-{mDark}" in k:
                        tempDict[k] = v
                        
            tempDict["QCD_Combined"] = sampleDict["QCD_Combined"]
    
            print("tempDict = ", tempDict)
            print()
            for region in ["HE", "HB"]:
    
                plot1DComparison(tempDict, f"leadJet_{region}_Pt", "Lead Jet p_{T} [GeV]" + f" ({region})", 0, 800,
                                 f"{args.outdir}/leadjetpt_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True)
                plot1DComparison(tempDict, f"leadJet_{region}_Eta", "Lead Jet #eta" + f" ({region})", -3, 3,
                                 f"{args.outdir}/leadjeteta_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}")
                plot1DComparison(tempDict, f"leadJet_{region}_Phi", "Lead Jet #phi" + f" ({region})", -3, 3,
                                 f"{args.outdir}/leadjetphi_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}")
                plot1DComparison(tempDict, f"leadJet_{region}_E", "Lead Jet E [GeV]" + f" ({region})", 0, 1000,
                                 f"{args.outdir}/leadjete_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True)
                plot1DComparison(tempDict, f"leadJet_{region}_nConstituent", "Lead Jet Number of Constituents" + f" ({region})", 0, 50,
                                 f"{args.outdir}/leadjetnconstituent_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}")
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_frac0", "Lead Jet Fraction of RecHits with auxTDC=0" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCfrac0_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_frac1", "Lead Jet Fraction of RecHits with auxTDC=1" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCfrac1_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_frac2", "Lead Jet Fraction of RecHits with auxTDC=2" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCfrac2_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_frac3", "Lead Jet Fraction of RecHits with auxTDC=3" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCfrac3_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_energyWeightedFrac0", "Lead Jet Energy-Weighted Fraction of RecHits with auxTDC=0" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCEnergyWeightedFrac0_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_energyWeightedFrac1", "Lead Jet Energy-Weighted Fraction of RecHits with auxTDC=1" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCEnergyWeightedFrac1_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_energyWeightedFrac2", "Lead Jet Energy-Weighted Fraction of RecHits with auxTDC=2" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCEnergyWeightedFrac2_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_auxTDC3_energyWeightedFrac3", "Lead Jet Energy-Weighted Fraction of RecHits with auxTDC=3" + f" ({region})", 0, 1,
                                 f"{args.outdir}/leadjet_auxTDCEnergyWeightedFrac3_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True, normStart=1)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_energy_avg", "Lead Jet Average RecHit Energy [GeV]" + f" ({region})", 0, 1, 
                                 f"{args.outdir}/leadjet_energyavg_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_energy_max", "Lead Jet Maximum RecHit Energy [GeV]" + f" ({region})", 0, 10, 
                                 f"{args.outdir}/leadjet_energymax_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True)
                plot1DComparison(tempDict, f"leadJet_{region}_JetRecHit_energy_min", "Lead Jet Minimum RecHit Energy [GeV]" + f" ({region})", 0, 0.5, 
                                 f"{args.outdir}/leadjet_energymin_{region}_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True)
                plot1DComparison(tempDict, f"recHit_{region}_E_max", r"Maximum per-Event RecHit Energy [GeV]" + f" ({region})", 0, 10,
                                 f"{args.outdir}/recHit_{region}_E_max_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True)
                plot1DComparison(tempDict, f"recHit_{region}_E_min", r"Minimum per-Event RecHit Energy [GeV]" + f" ({region})", 0, 0.5,
                                 f"{args.outdir}/recHit_{region}_E_min_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True)
                plot1DComparison(tempDict, f"recHit_{region}_E_median", r"Median per-Event RecHit Energy [GeV]" + f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_E_median_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True)
                plot1DComparison(tempDict, f"recHit_{region}_E_mean", r"Mean per-Event RecHit Energy [GeV]" + f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_E_mean_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_frac0", r"Fraction of RecHit auxTDC 3 equal to 0" + f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_auxTDC3_frac0_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_frac1", r"Fraction of RecHit auxTDC 3 equal to 1" + f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_auxTDC3_frac1_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_frac2", r"Fraction of RecHit auxTDC 3 equal to 2" + f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_auxTDC3_frac2_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_frac3", r"Fraction of RecHit auxTDC 3 equal to 3" + f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_auxTDC3_frac3_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_energyWeightedFrac0", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 0" + f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_auxTDC3_energyWeightedFrac0_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_energyWeightedFrac1", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 1" + f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_auxTDC3_energyWeightedFrac1_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_energyWeightedFrac2", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 2" + f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_auxTDC3_energyWeightedFrac2_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)
                plot1DComparison(tempDict, f"recHit_{region}_auxTDC3_energyWeightedFrac3", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 3" + f" ({region})", 0, 1,
                                 f"{args.outdir}/recHit_{region}_auxTDC3_energyWeightedFrac3_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)
            plot1DComparison(tempDict, "leadJet_Pt", "Lead Jet p_{T} [GeV]", 0, 800, f"{args.outdir}/leadjetpt_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True)
            plot1DComparison(tempDict, "leadJet_Eta", "Lead Jet #eta", -3, 3, f"{args.outdir}/leadjeteta_mDark-{mDark}_mMed-{mMed}_{todaysDate}")
            plot1DComparison(tempDict, "leadJet_Phi", "Lead Jet #phi", -3, 3, f"{args.outdir}/leadjetphi_mDark-{mDark}_mMed-{mMed}_{todaysDate}")
            plot1DComparison(tempDict, "leadJet_E", "Lead Jet E [GeV]", 0, 1000, f"{args.outdir}/leadjete_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True)
            plot1DComparison(tempDict, "leadJet_nConstituent", "Lead Jet Number of Constituents", 0, 50, f"{args.outdir}/leadjetnconstituent_mDark-{mDark}_mMed-{mMed}_{todaysDate}")
            plot1DComparison(tempDict, f"leadJet_JetRecHit_auxTDC3_frac0", "Lead Jet Fraction of RecHits with auxTDC=0", 0, 1, f"{args.outdir}/leadjet_auxTDCfrac0_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True, normStart=1)
            plot1DComparison(tempDict, f"leadJet_JetRecHit_auxTDC3_frac1", "Lead Jet Fraction of RecHits with auxTDC=1", 0, 1, f"{args.outdir}/leadjet_auxTDCfrac1_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True, normStart=1)
            plot1DComparison(tempDict, f"leadJet_JetRecHit_auxTDC3_frac2", "Lead Jet Fraction of RecHits with auxTDC=2", 0, 1, f"{args.outdir}/leadjet_auxTDCfrac2_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True, normStart=1)
            plot1DComparison(tempDict, f"leadJet_JetRecHit_auxTDC3_frac3", "Lead Jet Fraction of RecHits with auxTDC=3", 0, 1, f"{args.outdir}/leadjet_auxTDCfrac3_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True, normStart=1)
            plot1DComparison(tempDict, f"leadJet_JetRecHit_energy_avg", "Lead Jet Average RecHit Energy [GeV]", 0, 1, f"{args.outdir}/leadjet_energyavg_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True)
            plot1DComparison(tempDict, f"leadJet_JetRecHit_energy_max", "Lead Jet Maximum RecHit Energy [GeV]", 0, 10, f"{args.outdir}/leadjet_energymax_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True)
            plot1DComparison(tempDict, f"leadJet_JetRecHit_energy_min", "Lead Jet Minimum RecHit Energy [GeV]", 0, 0.5, f"{args.outdir}/leadjet_energymin_mDark-{mDark}_mMed-{mMed}_{todaysDate}", log=True)
            plot1DComparison(tempDict, f"recHit_E_max", r"Maximum per-Event RecHit Energy [GeV]", 0, 10, f"{args.outdir}/recHit_E_max_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True)
            plot1DComparison(tempDict, f"recHit_E_min", r"Minimum per-Event RecHit Energy [GeV]", 0, 0.5, f"{args.outdir}/recHit_E_min_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True)
            plot1DComparison(tempDict, f"recHit_E_median", r"Median per-Event RecHit Energy [GeV]", 0, 1,f"{args.outdir}/recHit_E_median_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True)
            plot1DComparison(tempDict, f"recHit_E_mean", r"Mean per-Event RecHit Energy [GeV]", 0, 1, f"{args.outdir}/recHit_E_mean_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True)
            plot1DComparison(tempDict, f"recHit_auxTDC3_frac0", r"Fraction of RecHit auxTDC 3 equal to 0", 0, 1, f"{args.outdir}/recHit_auxTDC3_frac0_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)
            plot1DComparison(tempDict, f"recHit_auxTDC3_frac1", r"Fraction of RecHit auxTDC 3 equal to 1", 0, 1, f"{args.outdir}/recHit_auxTDC3_frac1_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)
            plot1DComparison(tempDict, f"recHit_auxTDC3_frac2", r"Fraction of RecHit auxTDC 3 equal to 2", 0, 1, f"{args.outdir}/recHit_auxTDC3_frac2_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)
            plot1DComparison(tempDict, f"recHit_auxTDC3_frac3", r"Fraction of RecHit auxTDC 3 equal to 3", 0, 1, f"{args.outdir}/recHit_auxTDC3_frac3_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)
            plot1DComparison(tempDict, f"recHit_auxTDC3_energyWeightedFrac0", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 0", 0, 1, f"{args.outdir}/recHit_auxTDC3_energyWeightedFrac0_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)
            plot1DComparison(tempDict, f"recHit_auxTDC3_energyWeightedFrac1", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 1", 0, 1, f"{args.outdir}/recHit_auxTDC3_energyWeightedFrac1_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)
            plot1DComparison(tempDict, f"recHit_auxTDC3_energyWeightedFrac2", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 2", 0, 1, f"{args.outdir}/recHit_auxTDC3_energyWeightedFrac2_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)
            plot1DComparison(tempDict, f"recHit_auxTDC3_energyWeightedFrac3", r"Energy-Weighted Fraction of RecHit auxTDC 3 equal to 3", 0, 1, f"{args.outdir}/recHit_auxTDC3_energyWeightedFrac3_mDark-{mDark}_mMed-{mMed}_{todaysDate}",log=True, normStart=1)

###################################################################################################
# RUN SCRIPT

if __name__ == "__main__":
    main()
