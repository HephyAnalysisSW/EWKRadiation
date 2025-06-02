#!/usr/bin/"env python
''' simple analysis script 
'''
#
# Standard imports and batch mode
#
#import ROOT
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
vector.register_awkward()
import warnings
warnings.filterwarnings("ignore", category=UserWarning)



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
    gROOT.SetBatch(True)
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



    

# Open ROOT file and get the Events tree
file_path = "/eos/vbc/experiments/cms/store/data/Run2018D/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v2/2430000/04942A65-FE75-0743-B0F9-87E3E249D7C3.root"
with uproot.open(file_path) as file:
    tree = file["Events"]

    #tkeys = tree.keys()
    #for tk in tkeys:
    #    if "Muon" in tk:
    #        print(tk)
    #exit()

    # Define branches to read (extend this list later as needed)
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
    data = tree.arrays(branches, library="ak")
    
    

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
    print(t_vec[0])
    print(t_vec[0][1]["pt"])
    print(t_vec[0][0].pt)
    

    print("How many muons?")
    print(ak.num(t_vec,axis=1)) #axis gibt dim der Verschachtelung an

    #Kombiniere all muonen innerhalb events
    pair_combinations = ak.combinations(t_vec, n_muons, axis=1)
    idx_pair_combinations = ak.argcombinations(t_vec, n_muons, axis=1)
    
 
    print(f"How many combinations of {n_muons}?")
    print(ak.num(pair_combinations,axis=1))
    
    print(muon_charges)
    
    #Total charge == 0
    charge_combinations = ak.combinations(muon_charges, n_muons, axis=1)
    charges_sum = sum(charge_combinations[str(i)] for i in range(n_muons))#aksum funktioniert nicht??
    print("charge comb:", charge_combinations)

    mask = charges_sum == 0
    print("charge comb:", charge_combinations)
    print("charges sum:", charges_sum)

    valid_pairs = pair_combinations[mask]
    idx_valid_pairs = idx_pair_combinations[mask]
    print(valid_pairs)
    

    start_time = time.time()
    inv_mass = ak.sum((valid_pairs[str(i)] for i in range(n_muons))).mass

    print(inv_mass)




    print("Invariant mass of valid pairs:", inv_mass)
    end_time = time.time()

    print("Invariant mass valid pairs:", inv_mass_valid)
    print(f"Time to compute invariant mass valid: {end_time - start_time:.4f} seconds")

    
    mass_all = ak.flatten(inv_mass_valid) #h1   
    #best_pair = np.argmin( np.abs(inv_mass - 91) ) 
    print("inv mass valid")
    print(inv_mass_valid)
    print(ak.min(inv_mass_valid,axis=1))
    print("diff")
    print(ak.min(np.abs(inv_mass_valid - 91), axis=1)) #ziehe jedem paar 91 ab und nimm kleinste Abweichung
    min_inv_mass = ak.min(np.abs(inv_mass_valid - 91), axis=1, keepdims = True) #behalte nur die beste kombination für Z
    best_pair = ak.argmin(np.abs(inv_mass_valid - 91), axis=1, keepdims = True)
    print(f"  best_pair: ",best_pair)
    print(f"  Best invariant mass combination overall: ",inv_mass_valid[best_pair])
    ##print(f"  Pair of muons: ",idx_valid_pairs)
    ##print(f"  Best pair of muons: ",idx_valid_pairs[best_pair])
    mass_closeZ = ak.flatten(inv_mass_valid[best_pair])
    


    ##Apply a selection on the invariant mass and filter away the events very far from the Z peak
    inv_mass_sel = (inv_mass_valid > (91-50)) & (inv_mass_valid < (91+50))
    #print(inv_mass[inv_mass_sel])
    #exit()
    #idx_valid_pairs_sel = idx_valid_pairs[inv_mass_sel]
    #best_pair = np.argmin( np.abs(inv_mass[inv_mass_sel] - 91) )
    best_pair = ak.argmin(np.abs(inv_mass_valid[inv_mass_sel] - 91), axis=1, keepdims = True)
    print(f"  Best invariant mass combination in Z mass window: ",inv_mass_valid[inv_mass_sel])
    print(f"  Best invariant mass combination in Z mass window, 2nd method: ",inv_mass_valid[best_pair])
    print(f"  Best pair of muons in Z mass window: ",best_pair)
    #mass_window = np.array(inv_mass[inv_mass_sel][best_pair]).flatten()
    mass_window = ak.flatten(inv_mass_valid[inv_mass_sel])
    

    #print(sum(len(v) for event in selected_event for v in event.values()))
    #print(selected_event)
#def create_histogram(mass_values, hist_name, bins, x_min, x_max, file_name, hist_title=""):
    create_histogram(mass_all, f"h1{n_muons}", 100, 0, 200, f"{n_muons}combination_masses_all")
    if n_muons == 2:
        create_histogram(mass_closeZ, f"h2{n_muons}", 100, 0, 200, f"{n_muons}combination_mass_closeZ")
        create_histogram(mass_window, f" h3{n_muons}")

    print(f"{n_muons}combinationFilled histo")
    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)
    
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
