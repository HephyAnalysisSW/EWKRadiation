# pip install uproot numpy matplotlib
import uproot, numpy as np, matplotlib.pyplot as plt
from pathlib import Path

FILES = [
    ("AllRuns_6muon.root",         "6"),
    ("AllRuns_6muon_xmax500.root", "6_xmax500"),
]
TARGET_BINS = 15

def rebin_to_fixed_bins(values, edges, nbins):
    centers = 0.5*(edges[:-1] + edges[1:])
    data = np.repeat(centers, values.astype(int))
    xmin, xmax = edges[0], edges[-1]
    new_edges = np.linspace(xmin, xmax, nbins+1)
    counts, _ = np.histogram(data, bins=new_edges)
    return counts, new_edges

def rootish_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("bottom","left"):
        ax.spines[s].set_linewidth(1.4)
    ax.tick_params(axis="both", direction="in", length=7, width=1.4)
    ax.tick_params(axis="both", which="minor", direction="in", length=4, width=1.0)
    ax.minorticks_on()
    ax.grid(False)

def stats_box(ax, title, counts, edges):
    n = int(np.sum(counts))
    centers = 0.5*(edges[:-1] + edges[1:])
    # gewichtete Kennzahlen (wie bei ROOT-Hist)
    mean = np.average(centers, weights=counts) if n>0 else np.nan
    var  = np.average((centers-mean)**2, weights=counts) if n>0 else np.nan
    std  = np.sqrt(var)
    text = (f"{title}\n"
            f"Entries   {n}\n"
            f"Mean    {mean:6.1f}\n"
            f"Std Dev {std:6.1f}")
    ax.text(
        0.98, 0.98, text, transform=ax.transAxes,
        ha="right", va="top",
        fontsize=11,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="black", lw=1.0)
    )

for filename, hname in FILES:
    with uproot.open(filename) as f:
        h = f[hname]
        vals, edges = h.to_numpy()
        counts, new_edges = rebin_to_fixed_bins(vals, edges, TARGET_BINS)
        centers = 0.5*(new_edges[:-1] + new_edges[1:])
        yerr = np.sqrt(counts)

        fig, ax = plt.subplots(figsize=(8.4,5.6))
        rootish_axes(ax)

        # Outline (blau) – wie ROOT "HIST"
        ax.step(new_edges, np.r_[counts, counts[-1] if len(counts) else 0],
                where="post", lw=1.4, color="#2b4bf2")

        # Punkte + Fehlerbalken (schwarz)
        ax.errorbar(centers, counts, yerr=yerr, fmt="o",
                    ms=5.5, mfc="black", mec="black",
                    ecolor="black", elinewidth=1.2, capsize=0)

        ax.set_xlabel("Invariant mass (GeV)")
        ax.set_ylabel("Events")
        ax.set_title(hname, loc="right", fontsize=12, pad=8)
        stats_box(ax, hname, counts, new_edges)

        # y-Achse hübsch: ganzzahlige Ticks
        ymax = max( (counts.max() if counts.size else 1) * 1.25, 1 )
        ax.set_ylim(0, ymax)
        ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

        out_png = f"{Path(filename).stem}_{hname}_15bins_rootlook.png"
        out_pdf = f"{Path(filename).stem}_{hname}_15bins_rootlook.pdf"
        fig.tight_layout()
        fig.savefig(out_png, dpi=200)
        fig.savefig(out_pdf)   # scharf für Paper
        plt.close(fig)
        print("gespeichert:", out_png, "und", out_pdf)
