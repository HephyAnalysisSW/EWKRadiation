import ROOT
import glob
from array import array
import random

# Alle Root-Dateien durchgehen
for filename in glob.glob("AllRuns_*.root"):
    print(f"Öffne {filename} ...")
    f = ROOT.TFile.Open(filename)

    # Keys in der Datei durchgehen
    for key in f.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom("TH1"):  # falls es ein Histogramm ist
            xs, ys = [], []

            # Alle Bins durchgehen
            for i in range(1, obj.GetNbinsX()+1):
                cnt = int(obj.GetBinContent(i))
                xc  = obj.GetBinCenter(i)
                for _ in range(cnt):
                    xs.append(xc)
                    # optional leicht jitter in y, sonst liegen alle exakt auf 0
                    ys.append(0.0 + random.gauss(0,0.02))

            if xs:
                g = ROOT.TGraph(len(xs), array('d', xs), array('d', ys))
                g.SetMarkerStyle(20)   # volle Kreise
                g.SetMarkerSize(0.6)
                g.SetMarkerColor(ROOT.kBlack)

                c = ROOT.TCanvas()
                g.Draw("AP")
                g.GetYaxis().SetRangeUser(-0.5,0.5)  # Achse "flach"
                g.GetYaxis().SetNdivisions(0)        # y-Achse ohne ticks
                g.SetTitle(f"{obj.GetName()};Wert;")

                outname = filename.replace(".root", f"_{obj.GetName()}_rug.png")
                c.SaveAs(outname)
                print(f"  -> gespeichert: {outname}")

        elif obj.InheritsFrom("TCanvas"):  # falls schon ein Canvas gespeichert ist
            outname = filename.replace(".root", f"_{obj.GetName()}.png")
            obj.SaveAs(outname)
            print(f"  -> gespeichert: {outname}")

    f.Close()
