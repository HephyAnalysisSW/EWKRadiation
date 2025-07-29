#!/usr/bin/"env python
''' simple analysis script 
'''
#
# Standard imports and batch mode
#
#import ROOT
import subprocess
import time
from ROOT import gROOT, TH1, TH1D, TCanvas, TLegend, gStyle, TFile
from math import sqrt, cos, sin, pi, cosh
#from RootTools.core.standard import *

# Maybe we do this later 
#from EWKRadiation.tools.user import plot_directory
#from EWKRadiation.tools.objectSelection import getGoodMuons, alwaysTrue

#import argparse
#argParser = argparse.ArgumentParser(description = "Argument parser")
#argParser.add_argument('--logLevel',       action='store',      default='INFO',      nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
#argParser.add_argument('--plot_directory', action='store',      default='FourMuonInvariantMass')
#argParser.add_argument('--small',       action='store_true',                                                                        help="Run the file on a small sample (for test purpose), bool flag set to True if used" )
#args = argParser.parse_args()

import uproot
import awkward as ak
import numpy as np
import vector
import os
vector.register_awkward()
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

gROOT.SetBatch(True)
gStyle.SetOptStat(1111)

def create_histogram(mass_values, hist_name, bins, x_min, x_max, file_name, hist_title=""):
    hist = TH1D(hist_name, hist_title, bins, x_min, x_max)
    if len(mass_values) >0:
        for m in mass_values:
            hist.Fill(m)
    
    outfile = TFile(f"{file_name}.root","Recreate")
    outfile.cd()
    hist.Write(hist_name)
    outfile.Close()

    # create and safe canvas
    

    c = TCanvas(f"c_{file_name}", hist_title, 900,675)
    c.cd()
    hist.GetXaxis().SetTitle("Invariant mass (GeV)")
    hist.GetYaxis().SetTitle("Events")
    hist.SetMarkerStyle(20)
    hist.Draw("PE")
    c.Print(f"{file_name}.pdf")
    c.Print(f"{file_name}.png")
    c.Close()

    return hist


def merge_histograms_and_plot(n_muons):
    pattern = f"hist_{n_muons}muon_batch"
    merged_file = f"hist_{n_muons}muon_combined.root"
    hist_name = f"h1_{n_muons}"

    # Filter passende .root-Dateien
    input_files = [f for f in os.listdir(".") if f.startswith(pattern) and f.endswith(".root")]
    if not input_files:
        print(f"No input files for {n_muons} muons.")
        return

    # Merge mit hadd
    print(f"Merging: {input_files}")
    subprocess.run(["hadd", "-f", merged_file] + input_files)

    # Öffne gemergte Datei und zeichne Histogramm
    file = TFile.Open(merged_file)
    hist = file.Get(hist_name)

    if hist:
        c = TCanvas(f"c_combined_{n_muons}", f"{n_muons} Muons Combined", 900, 675)
        c.cd()
        hist.SetMarkerStyle(20)
        hist.GetXaxis().SetTitle("Invariant mass (GeV)")
        hist.GetYaxis().SetTitle("Events")
        hist.Draw("PE")
        c.Print(f"hist_{n_muons}muon_combined.pdf")
        c.Print(f"hist_{n_muons}muon_combined.png")
        c.Close()
        print(f"Saved merged plots for {n_muons} muons.")
    else:
        print(f"Histogram {hist_name} not found in {merged_file}.")

    file.Close()



def process_files(file_paths):
    print(f"Processing batch of files: {file_paths}")
    file_name_base = "batch" + str(abs(hash(",".join(file_paths)))) #basename gibt letzten /part
    
    
    branches = [
        "Muon_pt",
        "Muon_eta",
        "Muon_phi",
        "Muon_charge",
        "Muon_mediumPromptId", #weniger fake myonen z.B aus Hadron zerfällen (gibt Boolean zurück)
        "Muon_tightId", #noch strengere selection für muonen -> higgs 
        "Muon_pfRelIso04_all" #not considering muon near jets
        ]

    # Read branches
    data = uproot.concatenate({f: "Events" for f in file_paths}, expressions=branches, library="ak")
    print(f"Loaded {len(data)})events from {len(file_paths)}files")
    

#file_path = "/eos/vbc/experiments/cms/store/data/Run2018D/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v2/2430000/04942A65-FE75-0743-B0F9-87E3E249D7C3.root"



    
    

    for n_muons in [2,4,6]:
        # Apply selection criteria (string-based)
        #selection = ak.sum((data["Muon_pt"] > 5) & (data["Muon_mediumPromptId"] == 1), axis=1) >=2 # this should be >=2
        selection = ak.sum((data["Muon_pt"] > 5) & (data["Muon_tightId"] == 1) & (data["Muon_pfRelIso04_all"]<0.4) , axis=1)  >= n_muons
        selected_data = data[selection]
    
        #obj_selection = ((selected_data["Muon_pt"] > 5) & (selected_data["Muon_mediumPromptId"] == 1) )
        obj_selection = ((selected_data["Muon_pt"] > 5) & (selected_data["Muon_tightId"] == 1) & (selected_data["Muon_pfRelIso04_all"]<0.4))
    
    #build t_vec already, based on selected_data, that filters away events with less than 2 good muons
    #Liste von Listen. Einträge gefüllt mit dictionary 
        t_vec = ak.zip({
            "pt" : selected_data["Muon_pt"][obj_selection],
            "phi" : selected_data["Muon_phi"][obj_selection],
            "eta" : selected_data["Muon_eta"][obj_selection],
            "mass" : selected_data["Muon_charge"][obj_selection]*0,
        }, with_name="Momentum4D") 

        muon_charges = selected_data["Muon_charge"][obj_selection]

        print("Total number of events in file:", len(data))
        print("Number of selected events in file:", len(selected_data))
        # if len(selected_data) != 0:
        #     print(t_vec[0])
        #     print(t_vec[0][1]["pt"])
        #     print(t_vec[0][0].pt)
        

        print("How many muons?")
        print(ak.num(t_vec,axis=1)) #axis gibt dim der Verschachtelung an

        #Kombiniere all muonen innerhalb events
        pair_combinations = ak.combinations(t_vec, n_muons, axis=1)
        idx_pair_combinations = ak.argcombinations(t_vec, n_muons, axis=1)
        print(pair_combinations)
    
        print(f"How many combinations of {n_muons}?")
        print(ak.num(pair_combinations,axis=1))
        
        print(muon_charges)
        
        #Total charge == 0
        charge_combinations = ak.combinations(muon_charges, n_muons, axis=1)
        charges_sum = sum(charge_combinations[str(i)] for i in range(n_muons))

        mask = charges_sum == 0
        print("charge comb:", charge_combinations)
        print("charges sum:", charges_sum)

        valid_pairs = pair_combinations[mask]
        idx_valid_pairs = idx_pair_combinations[mask]
        #valid pairs already has empty arrays if event has not a single combination
        

        start_time = time.time()
        #inv_mass = ak.sum((valid_pairs[str(i)] for i in range(n_muons))).mass 
        #print(inv_mass)

        zero_momentum = ak.zeros_like(valid_pairs[str(0)])
        inv_mass = ( sum((valid_pairs[str(i)] for i in range(n_muons)), zero_momentum) ).mass
        print("sum well coded:")
        


        print("Invariant mass of valid pairs:", inv_mass)
        end_time = time.time()

        print(f"Time to compute invariant mass valid: {end_time - start_time:.4f} seconds")

        
        mass_all = ak.flatten(inv_mass) #h1  
        print("Number of entries:",len(mass_all)) 
        #best_pair = np.argmin( np.abs(inv_mass - 91) )
        
        out_file_base = f"hist_{n_muons}muon_batch{i}"  # i ist der Batch-Index
        create_histogram(mass_all, f"h1_{n_muons}", 100, 0, 200, out_file_base)
        print(f"Saved histogram for {n_muons} muons: {out_file_base}")
        
        
    
#     print(ak.min(inv_mass,axis=1))
#     print("diff")
#     print(ak.min(np.abs(inv_mass - 91), axis=1)) #ziehe jedem paar 91 ab und nimm kleinste Abweichung
#     min_inv_mass = ak.min(np.abs(inv_mass - 91), axis=1, keepdims = True) #behalte nur die beste kombination für Z
#     best_pair = ak.argmin(np.abs(inv_mass - 91), axis=1, keepdims = True)
#     print(f"  best_pair: ",best_pair)
#     print(f"  Best invariant mass combination overall: ",inv_mass[best_pair])
#     ##print(f"  Pair of muons: ",idx_valid_pairs)
#     ##print(f"  Best pair of muons: ",idx_valid_pairs[best_pair])
#     mass_closeZ = ak.flatten(inv_mass[best_pair])
    


#     ##Apply a selection on the invariant mass and filter away the events very far from the Z peak
#     inv_mass_sel = (inv_mass > (91-50)) & (inv_mass < (91+50))
#     #print(inv_mass[inv_mass_sel])
#     #exit()
#     #idx_valid_pairs_sel = idx_valid_pairs[inv_mass_sel]
#     #best_pair = np.argmin( np.abs(inv_mass[inv_mass_sel] - 91) )
#     best_pair = ak.argmin(np.abs(inv_mass[inv_mass_sel] - 91), axis=1, keepdims = True)
#     print(f"  Best invariant mass combination in Z mass window: ",inv_mass[inv_mass_sel])
#     print(f"  Best invariant mass combination in Z mass window, 2nd method: ",inv_mass[best_pair])
#     print(f"  Best pair of muons in Z mass window: ",best_pair)
#     #mass_window = np.array(inv_mass[inv_mass_sel][best_pair]).flatten()
#     mass_window = ak.flatten(inv_mass[inv_mass_sel])
    

#     #print(sum(len(v) for event in selected_event for v in event.values()))
#     #print(selected_event)
# #def create_histogram(mass_values, hist_name, bins, x_min, x_max, file_name, hist_title=""):
    
#     #if n_muons == 2:
#         #create_histogram(mass_closeZ, f"h2{n_muons}", 100, 0, 200, f"{n_muons}combination_mass_closeZ")
#         #create_histogram(mass_window, f" h3{n_muons}")

#     print(f"{n_muons}combinationFilled histo")
    
    #  #Draw a canvas
    # c1 = TCanvas("c1", "inv_mass_compare", 900, 675)
    # c1.cd()
    # h1.GetXaxis().SetTitle("Invariant mass (GeV)")
    # h1.GetYaxis().SetTitle("Events")
    # h1.SetLineColor(2)
    # h1.SetMarkerColor(2)
    # h1.SetMarkerStyle(20)
    # h1.Draw("PE")

    # h2.SetLineColor(4)
    # h2.SetMarkerColor(4)
    # h2.SetMarkerStyle(21)
    # h2.Draw("PE, sames")

    # h3.SetLineColor(8)
    # h3.SetMarkerColor(8)
    # h3.SetMarkerStyle(22)
    # h3.Draw("PE, sames")

    # legend = TLegend(0.6, 0.7, 0.88, 0.88)  # (x1, y1, x2, y2)
    # legend.AddEntry(h1, "No selection", "lep")
    # legend.AddEntry(h2, "Close to Z", "lep")
    # legend.AddEntry(h3, "Z window cut", "lep")
    # legend.Draw()

    # c1.Print("inv_mass_compare.pdf")
    # c1.Print("inv_mass_compare.png")
    # c1.Close()
    



#te_histogram(mass, noselection,100, 0, 200, noselection)   def makeM4l(event, sample):
#
#    # Get muons
#    muons      = getGoodMuons(event, collVars = muVars, mu_selector = lambda m:m['pt']>10 and m['mediumPromptId'] and m["pfRelIso03_all"]<0.3 )
#
#    event.m4l = float('nan')
#    m4l2 = 0
#    if len( muons ) ==4:
#
#        # select 2 positive and 2 negative charges
#        pdgIds = [ p['pdgId'] for p in muons ]
#        if pdgIds.count(+13) == 2 and pdgIds.count(-13)==2:
#            for i in range(1,4):
#                for j in range( i ):
#                    #print "indices i=%i,j=%i: adding %d to m4l2 => yiels %d" % (i,j,m4l2summand,m4l2)
#                    m4l2 += 27. #FIXME: Rosmarie. 
#
#    event.m4l = sqrt( m4l2 )
#    print "%d = invariant mass of 4 muons (2 pos, 2 neg)" % event.m4l
#   

file_list = [
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/0505F2F0-6F6D-0240-946A-9C2B65BFFA02.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/098ABF75-8D93-1B45-84A1-981613A831E7.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/20324905-EE60-9445-83C3-3F18A6C3E750.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/26A32EF4-689F-B742-BCA5-60E60FAE1D95.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/31F4F5F6-70B9-B341-BAF0-04FB8B240192.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/3867A6BC-12A8-C445-AE4F-7E9ABE5DC0D5.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/61A7CB73-0E0D-0A4A-BE72-4787CA196E3E.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/74F60BED-114B-BC4E-9FDA-CC36D414FF1D.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/7956F3EF-D268-564F-95EC-99DA4D5E7ADF.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/7D32F35F-A811-D844-A712-DBB6810370C2.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/7EE5C527-B391-934D-B59E-9AC40C2BC488.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/838866CF-18AE-8A45-8780-76F2E359D4C0.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/8AE6877F-91E6-914D-B346-09514695F15F.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/A9522567-DBC2-474B-B5CD-3667BB1ED13C.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/B62FFC66-F06B-B148-8B2C-E8AC12703159.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/C21C3C23-5723-D545-A3B4-74B0D5524A6F.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/C429D86F-5681-784E-99C0-0641C0D05EB4.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/D2CAFA6D-D802-2F48-BDDE-0B48CB1E0466.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/DF04EA1A-CC93-6A44-89D2-BCCB80F34A4E.root",
    "/eos/vbc/experiments/cms/store/data/Run2018B/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/E4B29848-0954-CA48-9123-9966DDE62D9A.root"
]

batch_size = 4

for i in range(0,len(file_list), batch_size):
    batch = file_list[i:i + batch_size]
    process_files(batch)
    

# Am Ende: Merged histos + PDFs + PNGs
for n in [2, 4, 6]:
    merge_histograms_and_plot(n)
