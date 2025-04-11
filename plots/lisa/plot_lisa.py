#!/usr/bin/env python
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


#Define a ROOT histogram
h = TH1D('h','', 100,0,200)

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
    obj_selection = ((data["Muon_pt"] > 5) & (data["Muon_mediumPromptId"] == 1) )
    selection = ak.sum((data["Muon_pt"] > 5) & (data["Muon_mediumPromptId"] == 1), axis=1) >= 2
    
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

    # Loop over selected events and access muon information
    for i_event in range(len(data["Muon_pt"])):
        if not selection[i_event]:
            continue

        print(f"Event {i_event}:")

        muon_pts = data["Muon_pt"][obj_selection][i_event]
        muon_etas = data["Muon_eta"][obj_selection][i_event]
        muon_phis = data["Muon_phi"][obj_selection][i_event]
        muon_charges = data["Muon_charge"][obj_selection][i_event]
        muon_mediumPromptIds = data["Muon_mediumPromptId"][obj_selection][i_event]

        n_muons = len(muon_pts)

        ##Philipp's code:
        ## 
        if n_muons>1:
            muon_4vecs = [four_momentum(muon_pts[i], muon_etas[i], muon_phis[i]) for i in range(n_muons)]
            print(muon_4vecs)
            print(f"mproduct test: {minkowski_product(muon_4vecs[0],muon_4vecs[0])}")
            print(f"invmass test: {inv_mass_funct(muon_4vecs[0],muon_4vecs[1])}")
            for i in range(n_muons):
                for j in range(i + 1, n_muons):
                    m_inv = inv_mass_funct(muon_4vecs[i], muon_4vecs[j])
                    print(f"Invariant Mass of muon pair ({i},{j}): {m_inv}:")


        for i_muon in range(n_muons):
            print(f"  Muon {i_muon}: pt={muon_pts[i_muon]}, eta={muon_etas[i_muon]}, phi={muon_phis[i_muon]}, "
                  f"charge={muon_charges[i_muon]}, mediumPromptId={muon_mediumPromptIds[i_muon]}")

        ##Option 2:
        ##Form a Lorentz vector, one for each muon, but looking at one event at a time
        t_vec = ak.zip({
            "pt" : data["Muon_pt"][obj_selection][i_event],
            "phi" : data["Muon_phi"][obj_selection][i_event],
            "eta" : data["Muon_eta"][obj_selection][i_event],
            "mass" : data["Muon_charge"][obj_selection][i_event]*0,
        }, with_name="Momentum4D")

        #print(ak.num(t_vec,axis=0))#this time I need axis=0, since there is only one dimension (number of muons)
        pair_combinations = ak.combinations(t_vec, 2, axis=0)
        idx_pair_combinations = ak.argcombinations(t_vec, 2, axis=0)
        #print(pair_combinations)
        #print(idx_pair_combinations)
        #print(ak.num(pair_combinations,axis=1))

        #Print out only the first 5 events
        #for pair in pair_combinations[0:5]:
        #    l = ak.num(pair,axis=0)
        #    print("Formed pairs: ",l)
        #    for pl in range(l):
        #        print(pl,pair[pl], (pair[pl]['0'] + pair[pl]['1']).mass)

        inv_mass = ak.Array([(pair['0'] + pair['1']).mass for pair in pair_combinations])

        ##No selections on the invariant masses
        print(f"  Invariant mass combinations per event: ",inv_mass)
        print(f"  Formed pair of muons: ",idx_pair_combinations)
        mass = np.array(inv_mass).flatten()
        #print(mass)

        ##If I can form more than one combination, I want the pair giving the closest invariant mass to the Z boson (91 GeV)
        best_pair = np.argmin( np.abs(inv_mass - 91) )
        print(f"  Best invariant mass combination: ",inv_mass[best_pair])
        print(f"  Best pair of muons: ",idx_pair_combinations[best_pair])
        mass = np.array(inv_mass[best_pair]).flatten()
        #print(mass)

        ##Apply a selection on the invariant mass and filter away the events very far from the Z peak
        inv_mass_sel = (inv_mass > (91-50)) & (inv_mass < (91+50))
        idx_pair_combinations_sel = idx_pair_combinations[inv_mass_sel]
        best_pair = np.argmin( np.abs(inv_mass[inv_mass_sel] - 91) )
        print(f"  Best invariant mass combination: ",inv_mass[inv_mass_sel][best_pair])
        print(f"  Best pair of muons: ",idx_pair_combinations_sel[best_pair])
        mass = np.array(inv_mass[inv_mass_sel][best_pair]).flatten()
        #print(mass)

        #mass can be either an array or a number; put each element in the histogram
        if len(mass)>0:
            for m in mass:
                h.Fill( m ) 

        ##if (mass!=None) or (ak.num(mass,axis=0)!=0):

        # Example: break after 5 events for brevity
        if i_event >= 50:#00:
            break

    print("Filled histo")
    gROOT.SetBatch(True)
    #gStyle.SetOptStat(0)

    #Draw a canvas
    c1 = TCanvas("c1", "inv_mass", 900, 675)
    c1.cd()
    h.GetXaxis().SetTitle("Invariant mass (GeV)")
    h.GetYaxis().SetTitle("Events")
    h.SetMarkerStyle(20)
    h.Draw("PE")
    c1.Print("inv_mass.pdf")
    c1.Print("inv_mass.png")
    c1.Close()

    #Write histogram in a root file
    #Warning! use a different label to avoid overwriting
    label = "_test"
    outfile = TFile("inv_mass"+label+".root","RECREATE")
    outfile.cd()
    h.Write("h")
    outfile.Close()
    print("Written ","inv_mass"+label+".root")



#def makeM4l(event, sample):
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
