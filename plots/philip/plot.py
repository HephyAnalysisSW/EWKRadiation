#!/usr/bin/env python
''' simple analysis script 
'''
#
# Standard imports and batch mode
#
#import ROOT
#ROOT.gROOT.SetBatch(True)
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
    selection = ak.sum((data["Muon_pt"] > 5) & (data["Muon_mediumPromptId"] == 1), axis=1) >= 4

    # Loop over selected events and access muon information
    for i_event in range(len(data["Muon_pt"])):
        if not selection[i_event]:
            continue

        print(f"Event {i_event}:")

        muon_pts = data["Muon_pt"][i_event]
        muon_etas = data["Muon_eta"][i_event]
        muon_phis = data["Muon_phi"][i_event]
        muon_charges = data["Muon_charge"][i_event]
        muon_mediumPromptIds = data["Muon_mediumPromptId"][i_event]

        n_muons = len(muon_pts)

        for i_muon in range(n_muons):
            print(f"  Muon {i_muon}: pt={muon_pts[i_muon]}, eta={muon_etas[i_muon]}, phi={muon_phis[i_muon]}, "
                  f"charge={muon_charges[i_muon]}, mediumPromptId={muon_mediumPromptIds[i_muon]}")

        # Example: break after 5 events for brevity
        if i_event >= 4:
            break

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
