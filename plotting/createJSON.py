import json

fileDict = {
    "QCD_Bin-Pt-30to50": "./data/fileList_Run3EMJJetStudy_QCD_Bin-Pt-30to50.txt",
    "QCD_Bin-Pt-50to80": "./data/fileList_Run3EMJJetStudy_QCD_Bin-Pt-50to80.txt",
    "QCD_Bin-Pt-80to120": "./data/fileList_Run3EMJJetStudy_QCD_Bin-Pt-80to120.txt",
    "QCD_Bin-Pt-120to170": "./data/fileList_Run3EMJJetStudy_QCD_Bin-Pt-120to170.txt",
    "QCD_Bin-Pt-170to300": "./data/fileList_Run3EMJJetStudy_QCD_Bin-Pt-170to300.txt",
    "QCD_Bin-Pt-300to470": "./data/fileList_Run3EMJJetStudy_QCD_Bin-Pt-300to470.txt",       
    "QCD_Bin-Pt-470to600": "./data/fileList_Run3EMJJetStudy_QCD_Bin-Pt-470to600.txt",
    "QCD_Bin-Pt-600to800": "./data/fileList_Run3EMJJetStudy_QCD_Bin-Pt-600to800.txt",
    "EMJ_s-channel_mMed-2000_mDark-10_ctau-1000_unflavored-down": "./data/fileList_s-channel_mMed-2000_mDark-10_ctau-1000_unflavored-down.txt",
    "EMJ_s-channel_mMed-2000_mDark-10_ctau-100_unflavored-down": "./data/fileList_s-channel_mMed-2000_mDark-10_ctau-100_unflavored-down.txt",
    "EMJ_s-channel_mMed-2000_mDark-10_ctau-1_unflavored-down": "./data/fileList_s-channel_mMed-2000_mDark-10_ctau-1_unflavored-down.txt",
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-1000_unflavored-down": "./data/fileList_s-channel_mMed-1500_mDark-10_ctau-1000_unflavored-down.txt",
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-100_unflavored-down": "./data/fileList_s-channel_mMed-1500_mDark-10_ctau-100_unflavored-down.txt",
    "EMJ_s-channel_mMed-1500_mDark-10_ctau-1_unflavored-down": "./data/fileList_s-channel_mMed-1500_mDark-10_ctau-1_unflavored-down.txt",
    "EMJ_s-channel_mMed-1500_mDark-20_ctau-1000_unflavored-down": "./data/fileList_s-channel_mMed-1500_mDark-20_ctau-1000_unflavored-down.txt",
    "EMJ_s-channel_mMed-1500_mDark-20_ctau-100_unflavored-down": "./data/fileList_s-channel_mMed-1500_mDark-20_ctau-100_unflavored-down.txt",
    "EMJ_s-channel_mMed-1500_mDark-20_ctau-1_unflavored-down": "./data/fileList_s-channel_mMed-1500_mDark-20_ctau-1_unflavored-down.txt",
    "EMJ_s-channel_mMed-1000_mDark-10_ctau-1000_unflavored-down": "./data/fileList_s-channel_mMed-1000_mDark-10_ctau-1000_unflavored-down.txt",
    "EMJ_s-channel_mMed-1000_mDark-10_ctau-100_unflavored-down": "./data/fileList_s-channel_mMed-1000_mDark-10_ctau-100_unflavored-down.txt",
    "EMJ_s-channel_mMed-1000_mDark-10_ctau-1_unflavored-down": "./data/fileList_s-channel_mMed-1000_mDark-10_ctau-1_unflavored-down.txt",
    "EMJ_s-channel_mMed-1000_mDark-20_ctau-1000_unflavored-down": "./data/fileList_s-channel_mMed-1000_mDark-20_ctau-1000_unflavored-down.txt",
    "EMJ_s-channel_mMed-1000_mDark-20_ctau-100_unflavored-down": "./data/fileList_s-channel_mMed-1000_mDark-20_ctau-100_unflavored-down.txt",
    "EMJ_s-channel_mMed-1000_mDark-20_ctau-1_unflavored-down": "./data/fileList_s-channel_mMed-1000_mDark-20_ctau-1_unflavored-down.txt",
    "EMJ_s-channel_mMed-750_mDark-10_ctau-1000_unflavored-down": "./data/fileList_s-channel_mMed-750_mDark-10_ctau-1000_unflavored-down.txt",
    "EMJ_s-channel_mMed-750_mDark-10_ctau-100_unflavored-down": "./data/fileList_s-channel_mMed-750_mDark-10_ctau-100_unflavored-down.txt",
    "EMJ_s-channel_mMed-750_mDark-10_ctau-1_unflavored-down": "./data/fileList_s-channel_mMed-750_mDark-10_ctau-1_unflavored-down.txt",
    "EMJ_s-channel_mMed-750_mDark-20_ctau-1000_unflavored-down": "./data/fileList_s-channel_mMed-750_mDark-20_ctau-1000_unflavored-down.txt",
    "EMJ_s-channel_mMed-750_mDark-20_ctau-100_unflavored-down": "./data/fileList_s-channel_mMed-750_mDark-20_ctau-100_unflavored-down.txt",
    "EMJ_s-channel_mMed-750_mDark-20_ctau-1_unflavored-down": "./data/fileList_s-channel_mMed-750_mDark-20_ctau-1_unflavored-down.txt",
    "EMJ_s-channel_mMed-500_mDark-10_ctau-1000_unflavored-down": "./data/fileList_s-channel_mMed-500_mDark-10_ctau-1000_unflavored-down.txt",
    "EMJ_s-channel_mMed-500_mDark-10_ctau-100_unflavored-down": "./data/fileList_s-channel_mMed-500_mDark-10_ctau-100_unflavored-down.txt",
    "EMJ_s-channel_mMed-500_mDark-10_ctau-1_unflavored-down": "./data/fileList_s-channel_mMed-500_mDark-10_ctau-1_unflavored-down.txt",
    "EMJ_s-channel_mMed-500_mDark-20_ctau-1000_unflavored-down": "./data/fileList_s-channel_mMed-500_mDark-20_ctau-1000_unflavored-down.txt",
    "EMJ_s-channel_mMed-500_mDark-20_ctau-100_unflavored-down": "./data/fileList_s-channel_mMed-500_mDark-20_ctau-100_unflavored-down.txt",
    "EMJ_s-channel_mMed-500_mDark-20_ctau-1_unflavored-down": "./data/fileList_s-channel_mMed-500_mDark-20_ctau-1_unflavored-down.txt",
    "EMJ_s-channel_mMed-250_mDark-10_ctau-1000_unflavored-down": "./data/fileList_s-channel_mMed-250_mDark-10_ctau-1000_unflavored-down.txt",
    "EMJ_s-channel_mMed-250_mDark-10_ctau-100_unflavored-down": "./data/fileList_s-channel_mMed-250_mDark-10_ctau-100_unflavored-down.txt",
    "EMJ_s-channel_mMed-250_mDark-10_ctau-1_unflavored-down": "./data/fileList_s-channel_mMed-250_mDark-10_ctau-1_unflavored-down.txt",
    "EMJ_s-channel_mMed-250_mDark-20_ctau-1000_unflavored-down": "./data/fileList_s-channel_mMed-250_mDark-20_ctau-1000_unflavored-down.txt",
    "EMJ_s-channel_mMed-250_mDark-20_ctau-100_unflavored-down": "./data/fileList_s-channel_mMed-250_mDark-20_ctau-100_unflavored-down.txt",
    "EMJ_s-channel_mMed-250_mDark-20_ctau-1_unflavored-down": "./data/fileList_s-channel_mMed-250_mDark-20_ctau-1_unflavored-down.txt",
    "EMJ_s-channel_mMed-100_mDark-10_ctau-1000_unflavored-down": "./data/fileList_s-channel_mMed-100_mDark-10_ctau-1000_unflavored-down.txt",
    "EMJ_s-channel_mMed-100_mDark-10_ctau-100_unflavored-down": "./data/fileList_s-channel_mMed-100_mDark-10_ctau-100_unflavored-down.txt",
    "EMJ_s-channel_mMed-100_mDark-10_ctau-1_unflavored-down": "./data/fileList_s-channel_mMed-100_mDark-10_ctau-1_unflavored-down.txt",
    "EMJ_s-channel_mMed-100_mDark-20_ctau-1000_unflavored-down": "./data/fileList_s-channel_mMed-100_mDark-20_ctau-1000_unflavored-down.txt",
    "EMJ_s-channel_mMed-100_mDark-20_ctau-100_unflavored-down": "./data/fileList_s-channel_mMed-100_mDark-20_ctau-100_unflavored-down.txt",
    "EMJ_s-channel_mMed-100_mDark-20_ctau-1_unflavored-down": "./data/fileList_s-channel_mMed-100_mDark-20_ctau-1_unflavored-down.txt",
}

# Read the text files and add lines to json dict
jsonDict = {}
for key in fileDict.keys():
    with open(fileDict[key], "r", encoding="utf-8") as file:
        jsonDict[key] = [line.strip() for line in file.readlines()]  # Strip to remove newline characters

# Write to a JSON file
with open("filePaths.json", "w", encoding="utf-8") as json_file:
    json.dump(jsonDict, json_file, indent=4)

print("JSON file created successfully: filePaths.json")
