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



REQS = {
    2: dict(pt_min=10.0, lead_pt_min=20.0, iso_branch="Muon_pfRelIso03_all", iso_max=0.35,
            use_sip3d=True,  sip3d_max=4.0),
    4: dict(pt_min=10.0, lead_pt_min=20.0, iso_branch="Muon_pfRelIso03_all", iso_max=0.35,
            use_sip3d=True,  sip3d_max=4.0),
    6: dict(pt_min=5.0,  lead_pt_min=None,  iso_branch="Muon_pfRelIso04_all", iso_max=0.35,
            use_sip3d=False, sip3d_max=None),
}



gROOT.SetBatch(True)
gStyle.SetOptStat(1111)

def create_histogram(mass_values, hist_name, bins, x_min, x_max, file_name_base, suffix=""):
    hist = TH1D(hist_name, "", bins, x_min, x_max)
    if len(mass_values) >0:
        for m in mass_values:
            hist.Fill(m)
    
    outname = f"{file_name_base}{suffix}"
    outfile = TFile(f"{outname}.root","Recreate")
    outfile.cd()
    hist.Write(hist_name)
    outfile.Close()

    # create and safe canvas
    

    c = TCanvas(f"c_{outname}", "", 900,675)
    c.cd()
    hist.GetXaxis().SetTitle("Invariant mass (GeV)")
    hist.GetYaxis().SetTitle("Events / {} GeV".format((x_max - x_min) // bins))
    hist.SetMarkerStyle(20)
    hist.Draw("PE")
    c.Print(f"{outname}.pdf")
    c.Print(f"{outname}.png")
    c.Close()

    return hist


def merge_histograms_and_plot(n_muons, list_name, bin_label=""):
    tag = list_name.replace("_file_list.txt", "")
    pattern = f"{tag}_{n_muons}muon_batch"
    merged_file = f"{tag}_{n_muons}muon_combined{bin_label}.root"
    hist_name = f"{n_muons}{bin_label}"

    #going through verzeichnis
    if bin_label == "":
        input_files = [f for f in os.listdir(".") if f.startswith(pattern) and not "_xmax500" in f and f.endswith(".root")]
    else:
        input_files = [f for f in os.listdir(".") if f.startswith(pattern) and f.endswith(f"{bin_label}.root")]

    if not input_files:
        print(f"No input files for {n_muons} muons with tag {tag}.")
        return

    subprocess.run(["hadd", "-f", merged_file] + input_files) #zusammenführung mit hadd

    file = TFile.Open(merged_file)
    hist = file.Get(hist_name)

    if hist:
        c = TCanvas(f"c_combined_{n_muons}{bin_label}_{tag}", "", 900, 675)
        c.cd()
        hist.SetMarkerStyle(20)
        hist.GetXaxis().SetTitle("Invariant mass (GeV)")
        hist.GetYaxis().SetTitle("Events / {} GeV".format(4 if "xmax500" in bin_label else 2))
        hist.Draw("PE")
        c.Print(f"{tag}_{n_muons}muon_combined{bin_label}.pdf")
        c.Print(f"{tag}_{n_muons}muon_combined{bin_label}.png")
        c.Close()
    else:
        print(f"Histogram {hist_name} not found in {merged_file}.")

    file.Close()

    #delete batchfiles
    for f in input_files:
        try:
            os.remove(f)
            for ext in [".pdf", ".png"]:
                side_file = f.replace(".root", ext)
                if os.path.exists(side_file):
                    os.remove(side_file)
        except Exception as e:
            print(f"Could not delete file {f}: {e}")


def process_files(file_paths, file_list_name, batch_index):

    tag = file_list_name.replace("_file_list.txt", "")
    file_name_base = f"{tag}"  # für alle Output-Namen

    branches = [
        "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge",
        "Muon_mediumPromptId", "Muon_looseId",
        "Muon_pfRelIso04_all", "Muon_pfRelIso03_all",
        #"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
        "Muon_sip3d",
    ]

#     with uproot.open(file_paths[0]) as test_file:
#         available_branches = test_file["Events"].keys()
#         print("\n📂 Checking available branches in file:", file_paths[0])
#         for b in branches:
#             if b not in available_branches:
#                 print(f"❌ Branch not found: {b}")
#             else:
#                 print(f"✅ Found: {b}")
                
    data = uproot.concatenate({f: "Events" for f in file_paths}, expressions=branches, library="ak")
#   hlt = data["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"] == True


    #Selections
    # ------------------ Selections (pro n_muons konfigurierbar) ------------------
    for n_muons in [2, 4, 6]:
        req = REQS[n_muons]
        iso_all = data[req["iso_branch"]]

        # --- Event-Ebene ---
        mu_mask_all = (data["Muon_pt"] > req["pt_min"]) \
                    & (data["Muon_looseId"] == 1) \
                    & (iso_all < req["iso_max"])
        if req["use_sip3d"]:
            mu_mask_all = mu_mask_all & (data["Muon_sip3d"] < req["sip3d_max"])

        sel_all = ak.sum(mu_mask_all, axis=1) >= n_muons

        # Lead-Requirement nur, wenn definiert
        if req["lead_pt_min"] is None:
            sel_lead = True
        else:
            mu_mask_lead = (data["Muon_pt"] > req["lead_pt_min"]) \
                        & (data["Muon_looseId"] == 1) \
                        & (iso_all < req["iso_max"])
            if req["use_sip3d"]:
                mu_mask_lead = mu_mask_lead & (data["Muon_sip3d"] < req["sip3d_max"])
            sel_lead = ak.any(mu_mask_lead, axis=1)

        combined_selection = sel_all & sel_lead
        selected_data = data[combined_selection]
        iso_sel = selected_data[req["iso_branch"]]

        # --- Objekt-Ebene für Kombinationen ---
        obj_sel = (selected_data["Muon_pt"] > req["pt_min"]) \
                & (selected_data["Muon_looseId"] == 1) \
                & (iso_sel < req["iso_max"])
        if req["use_sip3d"]:
            obj_sel = obj_sel & (selected_data["Muon_sip3d"] < req["sip3d_max"])

        t_vec = ak.zip({
            "pt":   selected_data["Muon_pt"][obj_sel],
            "phi":  selected_data["Muon_phi"][obj_sel],
            "eta":  selected_data["Muon_eta"][obj_sel],
            "mass": selected_data["Muon_charge"][obj_sel] * 0,
        }, with_name="Momentum4D")

        muon_charges = selected_data["Muon_charge"][obj_sel]

        comb   = ak.combinations(t_vec, n_muons, axis=1)
        charge = ak.combinations(muon_charges, n_muons, axis=1)
        charge_sum = sum(charge[str(i)] for i in range(n_muons))
        valid = comb[charge_sum == 0]

        zero = ak.zeros_like(valid["0"])
        inv_mass = (sum((valid[str(i)] for i in range(n_muons)), zero)).mass

        inv_mass_selected = ak.flatten(inv_mass[inv_mass > 4])

        hname     = f"{n_muons}"
        file_base = f"{file_name_base}_{n_muons}muon_batch{batch_index}"
        create_histogram(inv_mass_selected, hname, 100, 0, 200, file_base, suffix="")

        hname_ext = f"{n_muons}_xmax500"
        create_histogram(inv_mass_selected, hname_ext, 125, 0, 500, file_base, suffix="_xmax500")


        
        
    
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

def load_file_list(file_path):
    with open(file_path) as f:
        return [line.strip() for line in f if line.strip()]


file_list_dir = "cleaned_file_lists"
batch_size = 4
all_tags = []

# >>> Ab welchem Run starten?
START_FROM_TAG = None   # exakt wie der Tag aus dem Dateinamen (z.B. Run2016F)
START_BATCH    = 0            # optional: Startindex innerhalb dieses Runs (0, 4, 8, ... bei batch_size=4)

# --- Schleife über alle *_file_list_cleaned.txt im Ordner ---
for file_list_name in sorted(os.listdir(file_list_dir)):
    if not file_list_name.endswith("_file_list_cleaned.txt"):
        continue

    tag = file_list_name.replace("_file_list_cleaned.txt", "")

    # Vorherige Runs überspringen, bis wir bei START_FROM_TAG sind
    if START_FROM_TAG is not None and tag != START_FROM_TAG:
        print(f"⏭️  Skip {tag}")
        continue

    print(f"\n>>> Processing file list: {file_list_name}")
    all_tags.append(tag)

    file_list_path = os.path.join(file_list_dir, file_list_name)
    file_list = load_file_list(file_list_path)

    # Ab bestimmtem Index innerhalb des Runs starten (0, 4, 8, ... bei batch_size=4)
    for i in range(START_BATCH, len(file_list), batch_size):
        batch = file_list[i:i + batch_size]
        print(f"\n📂 Processing batch {i} for tag {tag} with files:")
        for f in batch:
            print(f"   - {f}")
        process_files(batch, tag, i)

    # Run-spezifisches Merging
    for n in [2, 4, 6]:
        merge_histograms_and_plot(n, tag, bin_label="")
        merge_histograms_and_plot(n, tag, bin_label="_xmax500")

    # Danach wieder alle folgenden Runs ganz normal mitnehmen
    START_FROM_TAG = None
    START_BATCH = 0


# --- End-Merge über alle Runs ---
print("\n✅ Alle Einzel-Runs verarbeitet. Starte nun mit dem finalen Zusammenführen der kombinierten Dateien...\n")
for n in [2, 4, 6]:
    for suffix in ["", "_xmax500"]:
        pattern = f"{n}muon_combined{suffix}.root"
        input_files = [f for f in os.listdir(".") if f.endswith(pattern)]
        output_file = f"AllRuns_{n}muon{suffix}.root"

        if input_files:
            print(f"📦 Merging {len(input_files)} Dateien zu {output_file}")
            subprocess.run(["hadd", "-f", output_file] + input_files)
