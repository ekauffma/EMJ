import json
import uproot

keys_of_interest = [
    'QCD_Bin-Pt-30to50', 
    'QCD_Bin-Pt-50to80', 
    'QCD_Bin-Pt-80to120', 
    'QCD_Bin-Pt-120to170', 
    'QCD_Bin-Pt-170to300', 
    'QCD_Bin-Pt-300to470', 
    'QCD_Bin-Pt-470to600', 
    'QCD_Bin-Pt-600to800'
]

with open('filePaths.json') as f:
    d = json.load(f)

for key in keys_of_interest:
    nEvents = 0
    for filePath in d[key]:
        f = uproot.open(filePath)
        nEvents+=f["Events"].num_entries
        f.close()
    print("Dataset = ", key, ", NEvents = ", nEvents)
