#!/usr/bin/"env python
''' simple analysis script 
'''
#
# Standard imports and batch mode
#
#import ROOT
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


def four_momentum(pt,eta,phi):
    p_0 = pt * np.cosh(eta)
    p_1 = pt * np.cos(phi)
    p_2 = pt * np.sin(phi)
    p_3 = pt * np.sinh(eta)
    return np.array([p_0,p_1,p_2,p_3])

def minkowski_product(p1,p2):
    return  -p1[0] * p2[0] + np.dot(p1[1:],p2[1:])

def inv_mass_funct(p1,p2):
    mass_squared = 2 * minkowski_product(p1,p2)
    return np.sqrt(np.abs(mass_squared))

def create_histogram(mass_values, hist_name, bins, x_min, x_max, label, file_name = "inv_mass",hist_title=""):
    hist = TH1D(hist_name, hist_title, bins, x_min, x_max)
    if len(mass_values) >0:
        for m in mass_values:
            hist.Fill(m)
    
    outfile = TFile(f"{file_name}{label}.root","Recreate")
    outfile.cd()
    hist.Write("h")
    outfile.Close()

    # create and safe canvas
    gROOT.SetBatch(True)
    c = TCanvas(f"c_{label}", hist_title, 900,675)
    c.cd()
    hist.GetXaxis().SetTitle("Invariant mass (Gev)")
    hist.GetYaxis().SetTitle("Events")
    hist.SetMarkerStyle(20)
    hist.Draw("PE")
    c.Print(f"{file_name}{label}.pdf")
    c.Print(f"{file_name}{label}.png")
    c.Close()

    return hist



# Open ROOT file and get the Events tree
file_path = "/eos/vbc/experiments/cms/store/data/Run2018D/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v2/2430000/04942A65-FE75-0743-B0F9-87E3E249D7C3.root"
with uproot.open(file_path) as file:
    tree = file["Events"]

    # Define branches to read (extend this list later as needed)
    branches = [
        "Muon_pt",
        "Muon_eta",
        "Muon_phi",
        "Muon_charge",
        "Muon_mediumPromptId"
    ]

    # Read branches
    data = tree.arrays(branches, library="ak")
    
    


    # Apply selection criteria (string-based)
    selection = ak.sum((data["Muon_pt"] > 5) & (data["Muon_mediumPromptId"] == 1), axis=1) >=2 # this should be >=2
    selected_data = data[selection]
    
    obj_selection = ((selected_data["Muon_pt"] > 5) & (selected_data["Muon_mediumPromptId"] == 1) )
    
    print("Total number of events in file:", len(data))
    print("Number of selected events in file:", len(selected_data))
    
    
    '''
    ##Option 1:
    ##Form a Lorentz vector, one for each muon in each event
    t_vec = ak.zip({
        "pt" : data["Muon_pt"][obj_selection][selection],
        "phi" : data["Muon_phi"][obj_selection][selection],
        "eta" : data["Muon_eta"][obj_selection][selection],
        "mass" : data["Muon_charge"][obj_selection][selection]*0,
    }, with_name="Momentum4D")

    ##Disadvantage: forming all of these combinations is a very slow process!
    ##As an example: combinations of muons present in the first 5 events
    print(ak.num(t_vec,axis=1))
    pair_combinations = ak.combinations(t_vec[0:5], 2, axis=1)
    print(ak.num(pair_combinations,axis=1))
    ##Print out the pairs formed in the first 5 events
    for pair in pair_combinations[0:5]:
        l = ak.num(pair,axis=0)
        print("Formed pairs: ",l)
        for pl in range(l):
            print(pl,pair[pl], (pair[pl]['0'] + pair[pl]['1']).mass)

    inv_mass = ak.from_iter([ak.Array([ (pair['0'] + pair['1']).mass for pair in pairs]) for pairs in pair_combinations])
    print("Invariant mass:", inv_mass)
    '''
    masses_all = []
    masses_closeZ = []
    masses_in_window = []

    selected_event = []
    sorted_out = []
    # Loop over selected events and access muon information
    for i_event in range(len(selected_data["Muon_pt"])):
        
        print(f"Event {i_event}:")

        muon_pts = selected_data["Muon_pt"][obj_selection][i_event]
        muon_etas = selected_data["Muon_eta"][obj_selection][i_event]
        muon_phis = selected_data["Muon_phi"][obj_selection][i_event]
        muon_charges = selected_data["Muon_charge"][obj_selection][i_event]
        muon_mediumPromptIds = selected_data["Muon_mediumPromptId"][obj_selection][i_event]

        n_muons = len(muon_pts)
       
 
        for i_muon in range(n_muons):
            print(f"  Muon {i_muon}: pt={muon_pts[i_muon]}, eta={muon_etas[i_muon]}, phi={muon_phis[i_muon]}, "
                  f"charge={muon_charges[i_muon]}, mediumPromptId={muon_mediumPromptIds[i_muon]}")

        ##Option 2:
        ##Form a Lorentz vector, one for each muon, but looking at one event at a time
        t_vec = ak.zip({
            "pt" : selected_data["Muon_pt"][obj_selection][i_event],
            "phi" : selected_data["Muon_phi"][obj_selection][i_event],
            "eta" : selected_data["Muon_eta"][obj_selection][i_event],
            "mass" : selected_data["Muon_charge"][obj_selection][i_event]*0,
        }, with_name="Momentum4D")
        
        print(t_vec)
        #print(ak.num(t_vec,axis=0))#this time I need axis=0, since there is only one dimension (number of muons)
        pair_combinations = ak.combinations(t_vec, 2, axis=0)
        idx_pair_combinations = ak.argcombinations(t_vec, 2, axis=0)

        #combinations total charge = 0
        charges_sum = muon_charges[idx_pair_combinations["0"]] + muon_charges[idx_pair_combinations["1"]] # + muon_charges[idx_pair_combinations["2"]] + muon_charges[idx_pair_combinations["3"]]
        mask = charges_sum == 0

        if ak.sum(mask) == 0:
            #sorted_out.append({i_event : valid_pairs})
            print("Total charge criteria not fullfilled")
            continue
        
        valid_pairs = pair_combinations[mask]
        idx_valid_pairs = idx_pair_combinations[mask]

        inv_mass = ak.Array([(pair['0'] + pair['1']).mass for pair in valid_pairs])
        
        selected_event.append({i_event : valid_pairs})
      
        
        

        ##No selections on the invariant masses
        print(f"  Invariant mass combinations per valid event: ",inv_mass)
        print(f"  Formed pair of valid muons: ",idx_valid_pairs)
        mass_all = np.array(inv_mass).flatten()
        masses_all.extend(mass_all) #h1   
        #print(mass)
        ##If I can form more than one combination, I want the pair giving the closest invariant mass to the Z boson (91 GeV)
        best_pair = np.argmin( np.abs(inv_mass - 91) )
        print(f"  Best invariant mass combination: ",inv_mass[best_pair])
        print(f"  Best pair of muons: ",idx_valid_pairs[best_pair])
        mass_closeZ = np.array(inv_mass[best_pair]).flatten()
        masses_closeZ.extend(mass_closeZ)
       #print(mass)

        ##Apply a selection on the invariant mass and filter away the events very far from the Z peak
        inv_mass_sel = (inv_mass > (91-50)) & (inv_mass < (91+50))
        idx_valid_pairs_sel = idx_valid_pairs[inv_mass_sel]
        best_pair = np.argmin( np.abs(inv_mass[inv_mass_sel] - 91) )
        print(f"  Best invariant mass combination: ",inv_mass[inv_mass_sel][best_pair])
        print(f"  Best pair of muons: ",idx_valid_pairs_sel[best_pair])
        mass_window = np.array(inv_mass[inv_mass_sel][best_pair]).flatten()
        masses_in_window.extend(mass_window)
       #print(mass)
 

        ##if (mass!=None) or (ak.num(mass,axis=0)!=0):

        # Example: break after 5 events for brevity
        if i_event >= 500:
            break
        
    print(sum(len(v) for event in selected_event for v in event.values()))
    print(selected_event)

    h1 = create_histogram(masses_all, "h1", 100, 0, 200, "masses_all")
    h2 = create_histogram(masses_closeZ, "h2", 100, 0, 200, "masses_closeZ")
    h3 = create_histogram(masses_in_window, "h3", 100, 0, 200, "masses_in_window")

    #print(f"Collected Mass for h1: {masses_all}")
    #print(f"Collected Mass for h2: {masses_closeZ}")
    #print(f"Collected Mass for h3: {masses_in_window}")
    #print("Filled histo")
    gROOT.SetBatch(True)
    #gStyle.SetOptStat(0)
    
    #Draw a canvas
    c1 = TCanvas("c1", "inv_mass_compare", 900, 675)
    c1.cd()
    h1.GetXaxis().SetTitle("Invariant mass (GeV)")
    h1.GetYaxis().SetTitle("Events")
    h1.SetLineColor(2)
    h1.SetMarkerStyle(20)
    h1.Draw("PE")

    h2.SetLineColor(4)
    h2.SetMarkerStyle(21)
    h2.Draw("PE, sames")

    h3.SetLineColor(8)
    h3.SetMarkerStyle(22)
    h3.Draw("PE, sames")

    legend = TLegend(0.6, 0.7, 0.88, 0.88)  # (x1, y1, x2, y2)
    legend.AddEntry(h1, "No selection", "lep")
    legend.AddEntry(h2, "Close to Z", "lep")
    legend.AddEntry(h3, "Z window cut", "lep")
    legend.Draw()

    c1.Print("inv_mass_compare.pdf")
    c1.Print("inv_mass_compare.png")
    c1.Close()
    



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
