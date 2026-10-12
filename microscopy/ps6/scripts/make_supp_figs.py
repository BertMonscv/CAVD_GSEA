# Publication copy (2026-10-12): analysis and plotting logic unchanged.
# Historical output figure labels are mapped to current manuscript figures in README.
#!/usr/bin/env python3
"""Supplementary Figures S11-S13 for the p-S6 (Ser235/236) immunofluorescence dataset.

S11  Focus-quality control and the prospectively defined sensitivity analysis
S12  Cell-density and segmentation covariate checks
S13  Complete field-level dataset (nothing aggregated away)

House style locked for this manuscript: 170 mm (6.69 in) canvas, 600 dpi,
Liberation Sans, Type-42 embedded fonts, no bbox_inches="tight" (it would trim
the canvas below 170 mm), left/bottom spines only, exact P values.
"""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator
from scipy import stats
from PIL import Image

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8.5, "axes.titlesize": 9.5,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.5,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
    "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    "xtick.minor.size": 1.7, "ytick.minor.size": 1.7,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.labelpad": 3.0, "lines.solid_capstyle": "round",
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})

# --------------------------------------------------------------------------- paths
# Resolved relative to this script, so the deposit runs anywhere after cloning:
#
#     <root>/data/percell.csv, perfield.csv, PS6_scan.csv
#     <root>/scripts/<this file>
#     <root>/figures/          <- output, created if absent
#
# If that layout is not found, the script falls back to looking beside itself and
# then in the current working directory. Both locations can be overridden:
#     PS6_DATA=/path/to/data  PS6_FIGURES=/path/to/out  python3 <this file>
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def _find_data():
    env = os.environ.get("PS6_DATA")
    if env:
        return Path(env)
    for cand in (ROOT / "data", HERE / "data", HERE, ROOT, Path.cwd()):
        if (cand / "percell.csv").exists():
            return cand
    sys.exit("percell.csv not found. Put the data files in <root>/data/ or set PS6_DATA.")


DATA = _find_data()
OUTDIR = Path(os.environ.get("PS6_FIGURES", ROOT / "figures"))
OUTDIR.mkdir(parents=True, exist_ok=True)
print("data   :", DATA)
print("figures:", OUTDIR)

cell = pd.read_csv(DATA / "percell.csv")
fld = pd.read_csv(DATA / "perfield.csv")
scan = pd.read_csv(DATA / "PS6_scan.csv")

for c in ["nCells", "net_field", "dapiNucEdge", "fieldMode", "fieldMax", "ringMedian_field"]:
    fld[c] = pd.to_numeric(fld[c], errors="coerce")
ok = fld[fld.status == "ok"].copy()

# perfield.csv's ringMedian_field / net_field are built from each cell's ring MEDIAN, while
# every statistic in the paper uses the ring MEAN. Recompute the per-field value from
# percell.csv so the field-level panels match Section 3.5 and Table S5.
ok = ok.join(cell.groupby(["file", "series"]).net.median().rename("net_field_true"),
             on=["file", "series"])
ok["rep"] = ok["rep"].astype(str)
cell["rep"] = cell["rep"].astype(str)

G = ["CTL", "PM", "LEP", "RAPA"]
ALL = G + ["Sac"]
LAB = {"CTL": "CTL", "PM": "PM", "LEP": "PM+Lep", "RAPA": "PM+Lep\n+Rapa", "Sac": "2°-only"}
FLAT = {k: v.replace("\n", "") for k, v in LAB.items()}
GC = {"CTL": "#55606B", "PM": "#2E6DB4", "LEP": "#C0392B", "RAPA": "#1E8A5A", "Sac": "#AEB4B9"}
RC = {"1": "#4E79A7", "2": "#F28E2B", "3": "#59A14F"}
THR = 0.25

ok = ok.merge(scan[["file", "series_idx", "C0_med"]],
              left_on=["file", "series"], right_on=["file", "series_idx"], how="left")


def covmed(d):
    f = (d.groupby(["group", "rep", "file", "series"], as_index=False).net.median()
         if "net" in d else None)
    return f


def stack_from_cells(sel_fields):
    """per-cell -> per-field median -> per-coverslip median, restricted to sel_fields."""
    key = set(zip(sel_fields.file, sel_fields.series))
    c = cell[cell.group.isin(G)].copy()
    c = c[[(f, s) in key for f, s in zip(c.file, c.series)]]
    fm = c.groupby(["group", "rep", "file", "series"], as_index=False).net.median()
    return fm.groupby(["group", "rep"], as_index=False).net.median().pivot(
        index="rep", columns="group", values="net")[G]


CONTR = [("CTL", "PM"), ("PM", "LEP"), ("LEP", "RAPA")]


def pvals(piv):
    return {c: stats.ttest_rel(piv[c[0]], piv[c[1]]).pvalue for c in CONTR}


def pf(p):
    return "P < 0.001" if p < 0.001 else "P = %.3f" % p


def plabel(ax, ch, x=-0.185, y=1.048):
    ax.text(x, y, ch.lower(), transform=ax.transAxes, fontsize=9.5,
            fontweight="bold", ha="left", va="bottom")


def tidy(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(which="both", pad=1.8)


def save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(OUTDIR / (stem + "." + ext), dpi=600, facecolor="white")
    plt.close(fig)
    im = Image.open(OUTDIR / (stem + ".png")).convert("RGB")
    im.save(OUTDIR / (stem + ".tif"), compression="tiff_lzw", dpi=(600, 600))
    print("%-28s %s px = %.1f mm wide" % (stem, im.size, im.size[0] / 600 * 25.4))


def jitter(n, centre, half, rng):
    return centre + rng.uniform(-half, half, n)


rng = np.random.default_rng(20260806)

# ══════════════════════════════════════════════════════════════════ S11
fig = plt.figure(figsize=(6.69, 5.90))
gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.32,
                      left=0.112, right=0.984, top=0.945, bottom=0.090)

# --- A: focus index by group
ax = fig.add_subplot(gs[0, 0])
for i, g in enumerate(ALL):
    v = ok[ok.group == g].dapiNucEdge.values
    ax.scatter(jitter(v.size, i, 0.20, rng), v, s=6, facecolor=GC[g], alpha=0.65,
               lw=0.2, edgecolor="white", zorder=3)
    m = np.median(v)
    ax.plot([i - 0.32, i + 0.32], [m, m], color="black", lw=1.3, zorder=5,
            solid_capstyle="butt")
    ax.text(i, 0.585, "%.3f" % m, ha="center", fontsize=7.2)
ax.axhline(THR, color="#C0392B", lw=0.8, ls=(0, (4, 2)), zorder=2)
ax.text(-0.52, THR + 0.008, "sensitivity-analysis threshold (0.25)", fontsize=6.9,
        ha="left", va="bottom", color="#C0392B")
kw = stats.kruskal(*[ok[ok.group == g].dapiNucEdge.dropna() for g in G]).pvalue
ax.text(0.98, 0.02, "Kruskal–Wallis across the four\ntreatment groups: P = %.4f" % kw,
        transform=ax.transAxes, fontsize=7.2, va="bottom", ha="right")
ax.set_xticks(range(5)); ax.set_xticklabels([LAB[g] for g in ALL])
ax.set_xlim(-0.6, 4.6); ax.set_ylim(0.10, 0.63)
ax.set_ylabel("Focus index (DAPI nuclear gradient)")
ax.yaxis.set_minor_locator(AutoMinorLocator(2)); plabel(ax, "A"); tidy(ax)

# --- B: focus index vs DAPI background (the residual confound, disclosed)
ax = fig.add_subplot(gs[0, 1])
for g in ALL:
    s = ok[ok.group == g]
    ax.scatter(s.C0_med, s.dapiNucEdge, s=7, facecolor=GC[g], alpha=0.7, lw=0.2,
               edgecolor="white", label=FLAT[g], zorder=3)
d = ok.dropna(subset=["C0_med", "dapiNucEdge"])
r, p = stats.pearsonr(d.C0_med, d.dapiNucEdge)
R_FOCUS, P_FOCUS = r, p
b1, b0 = np.polyfit(d.C0_med, d.dapiNucEdge, 1)
xs = np.linspace(d.C0_med.min(), d.C0_med.max(), 50)
ax.plot(xs, b0 + b1 * xs, color="black", lw=0.9, zorder=4)
ax.text(0.97, 0.95, "r = %+.3f, %s\n(n = %d fields)" % (r, pf(p), len(d)),
        transform=ax.transAxes, ha="right", va="top", fontsize=7.3)
ax.set_xlabel("Field median DAPI intensity (a.u.)")
ax.set_ylabel("Focus index")
ax.legend(loc="lower left", frameon=False, handletextpad=0.2, borderpad=0.05,
          labelspacing=0.2, markerscale=1.3)
ax.xaxis.set_minor_locator(AutoMinorLocator(2)); ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plabel(ax, "B", x=-0.165); tidy(ax)

# --- C: main vs sensitivity, coverslip medians
ax = fig.add_subplot(gs[1, 0])
piv_all = stack_from_cells(ok)
piv_sen = stack_from_cells(ok[ok.dapiNucEdge >= THR])
p_all, p_sen = pvals(piv_all), pvals(piv_sen)
for i, g in enumerate(G):
    for j, (piv, dx) in enumerate([(piv_all, -0.17), (piv_sen, 0.17)]):
        for r_ in "123":
            ax.scatter(i + dx, piv.loc[r_, g], s=26, marker="D", facecolor=RC[r_],
                       edgecolor="black", lw=0.5, zorder=5)
        m = piv[g].mean()
        ax.plot([i + dx - 0.13, i + dx + 0.13], [m, m], color="black", lw=1.2, zorder=6)
    for r_ in "123":
        ax.plot([i - 0.17, i + 0.17], [piv_all.loc[r_, g], piv_sen.loc[r_, g]],
                color="#9AA0A6", lw=0.5, zorder=2)
ax.set_xticks(range(4)); ax.set_xticklabels([LAB[g] for g in G])
ax.set_xlim(-0.55, 3.55)
ax.set_ylabel("Background-subtracted intensity\n(per-coverslip median, a.u.)")
ax.text(0.02, 0.97, "left of each pair: all fields\nright: focus index ≥ 0.25",
        transform=ax.transAxes, fontsize=7.2, va="top")
ax.yaxis.set_minor_locator(AutoMinorLocator(2)); plabel(ax, "C"); tidy(ax)

# --- D: P values across thresholds
ax = fig.add_subplot(gs[1, 1])
ths = [0.0, 0.25, 0.28, 0.30]
series = {c: [] for c in CONTR}
kept = []
for t in ths:
    sub = ok[ok.dapiNucEdge >= t] if t > 0 else ok
    kept.append(len(sub))
    pv = pvals(stack_from_cells(sub))
    for c in CONTR:
        series[c].append(pv[c])
STY = {("CTL", "PM"): ("#55606B", "o"), ("PM", "LEP"): ("#C0392B", "s"),
       ("LEP", "RAPA"): ("#1E8A5A", "^")}
for c in CONTR:
    col, mk = STY[c]
    ax.plot(range(len(ths)), series[c], color=col, lw=1.1, marker=mk, ms=4,
            mec="black", mew=0.4, label="%s vs %s" % (FLAT[c[0]], FLAT[c[1]]), zorder=4)
ax.axhline(0.05, color="#9AA0A6", lw=0.7, ls=":", zorder=2)
ax.text(-0.42, 0.056, "P = 0.05", fontsize=6.9, ha="left", va="bottom", color="#6b7075")
ax.set_yscale("log"); ax.set_ylim(2e-4, 1)
ax.set_xticks(range(len(ths)))
ax.set_xticklabels(["none\n(%d)" % kept[0]] + ["≥ %.2f\n(%d)" % (t, k)
                                               for t, k in zip(ths[1:], kept[1:])])
ax.set_xlabel("Focus-index threshold (fields retained)")
ax.set_ylabel("P value (paired t, Holm-uncorrected)")
ax.legend(loc="lower left", frameon=False, handlelength=1.4, handletextpad=0.4,
          borderpad=0.05, labelspacing=0.25)
plabel(ax, "D", x=-0.165); tidy(ax)
save(fig, "FigureS10_pS6_focusQC")

# ══════════════════════════════════════════════════════════════════ S12
fig = plt.figure(figsize=(6.69, 5.90))
gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.32,
                      left=0.112, right=0.984, top=0.945, bottom=0.090)

# --- A: net vs nuclei per field
ax = fig.add_subplot(gs[0, 0])
txt = []
for g in G:
    s = ok[ok.group == g]
    ax.scatter(s.nCells, s.net_field_true, s=8, facecolor=GC[g], alpha=0.7, lw=0.2,
               edgecolor="white", zorder=3)
    b1, b0 = np.polyfit(s.nCells, s.net_field_true, 1)
    xs = np.linspace(s.nCells.min(), s.nCells.max(), 20)
    ax.plot(xs, b0 + b1 * xs, color=GC[g], lw=0.9, ls=(0, (3, 2)), zorder=4)
    r, p = stats.pearsonr(s.nCells, s.net_field_true)
    txt.append("%-11s r = %+.2f, P = %.2f" % (FLAT[g], r, p))
for k, g in enumerate(G):
    ax.text(0.02, 0.975 - 0.062 * k, txt[k], transform=ax.transAxes, ha="left", va="top",
            fontsize=7.1, family="monospace", color=GC[g])
ax.set_ylim(120, 1900)
ax.set_xlabel("Nuclei per field"); ax.set_ylabel("Background-subtracted intensity\n(per-field median, a.u.)")
ax.xaxis.set_minor_locator(AutoMinorLocator(2)); ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plabel(ax, "A"); tidy(ax)

# --- B: field background vs nuclei
ax = fig.add_subplot(gs[0, 1])
txt = []
for g in G:
    s = ok[ok.group == g]
    ax.scatter(s.nCells, s.fieldMode, s=8, facecolor=GC[g], alpha=0.7, lw=0.2,
               edgecolor="white", zorder=3)
    b1, b0 = np.polyfit(s.nCells, s.fieldMode, 1)
    xs = np.linspace(s.nCells.min(), s.nCells.max(), 20)
    ax.plot(xs, b0 + b1 * xs, color=GC[g], lw=0.9, ls=(0, (3, 2)), zorder=4)
    r, p = stats.pearsonr(s.nCells, s.fieldMode)
    txt.append("%-11s r = %+.2f, P = %.3f" % (FLAT[g], r, p))
for k, g in enumerate(G):
    ax.text(0.02, 0.975 - 0.062 * k, txt[k], transform=ax.transAxes, ha="left", va="top",
            fontsize=7.1, family="monospace", color=GC[g])
ax.set_ylim(500, 2350)
ax.set_xlabel("Nuclei per field"); ax.set_ylabel("Field modal background (a.u.)")
ax.xaxis.set_minor_locator(AutoMinorLocator(2)); ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plabel(ax, "B", x=-0.165); tidy(ax)

# --- C: ring area per group
ax = fig.add_subplot(gs[1, 0])
for i, g in enumerate(ALL):
    v = cell[cell.group == g].ringArea_um2.values
    parts = ax.violinplot([v], positions=[i], widths=0.72, showextrema=False, showmedians=False)
    for b in parts["bodies"]:
        b.set_facecolor(GC[g]); b.set_alpha(0.30); b.set_edgecolor(GC[g]); b.set_linewidth(0.6)
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    ax.plot([i, i], [q1, q3], color="black", lw=2.4, solid_capstyle="butt", zorder=5)
    ax.plot([i], [med], marker="o", ms=3.2, mfc="white", mec="black", mew=0.7, zorder=6)
    ax.text(i, 640, "%.0f" % med, ha="center", fontsize=7.2)
ax.set_xticks(range(5)); ax.set_xticklabels([LAB[g] for g in ALL])
ax.set_xlim(-0.6, 4.6); ax.set_ylim(120, 690)
ax.set_ylabel("Perinuclear ring area per cell (µm²)")
ax.yaxis.set_minor_locator(AutoMinorLocator(2)); plabel(ax, "C"); tidy(ax)

# --- D: per-cell net vs ring area
ax = fig.add_subplot(gs[1, 1])
for g in ["Sac", "RAPA", "CTL", "PM", "LEP"]:
    s = cell[cell.group == g]
    ax.scatter(s.ringArea_um2, s.net, s=1.2, facecolor=GC[g], alpha=0.30, lw=0,
               rasterized=True, zorder=3)
r, p = stats.pearsonr(cell.ringArea_um2, cell.net)
ax.text(0.97, 0.97, "all cells pooled: r = %+.3f\n(n = %d)" % (r, len(cell)),
        transform=ax.transAxes, ha="right", va="top", fontsize=7.3)
ax.set_yscale("log"); ax.set_ylim(30, 4000)
ax.set_xlabel("Perinuclear ring area (µm²)")
ax.set_ylabel("Background-subtracted intensity per cell (a.u.)")
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
plabel(ax, "D", x=-0.165); tidy(ax)
save(fig, "FigureS11_pS6_covariates")

# ══════════════════════════════════════════════════════════════════ S13
order = [("CTL1.lif", "CTL"), ("CTL2.lif", "CTL"), ("CTL3.lif", "CTL"),
         ("PM1.lif", "PM"), ("PM2.lif", "PM"), ("PM3.lif", "PM"),
         ("LEP1.lif", "LEP"), ("LEP2.lif", "LEP"), ("LEP3.lif", "LEP"),
         ("RAPA1.lif", "RAPA"), ("RAPA2.lif", "RAPA"), ("RAPA3.lif", "RAPA"),
         ("PS6-Sac.lif", "Sac")]
names = [f.replace(".lif", "") for f, _ in order]

fig = plt.figure(figsize=(6.69, 5.90))
gs = fig.add_gridspec(2, 2, hspace=0.60, wspace=0.32,
                      left=0.105, right=0.982, top=0.940, bottom=0.135)


def bycov(ax, col, ylab, logy=False):
    for i, (f, g) in enumerate(order):
        v = ok[ok.file == f][col].values
        ax.scatter(jitter(v.size, i, 0.19, rng), v, s=7, facecolor=GC[g], alpha=0.75,
                   lw=0.2, edgecolor="white", zorder=3)
        m = np.median(v)
        ax.plot([i - 0.33, i + 0.33], [m, m], color="black", lw=1.2, zorder=5,
                solid_capstyle="butt")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(names, rotation=60, ha="right", fontsize=6.8)
    ax.set_xlim(-0.7, len(order) - 0.3)
    ax.set_ylabel(ylab)
    if logy:
        ax.set_yscale("log")
    else:
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    tidy(ax)


ax = fig.add_subplot(gs[0, 0]); bycov(ax, "net_field_true",
    "Background-subtracted intensity\n(per-field median, a.u.)"); plabel(ax, "A")
ax = fig.add_subplot(gs[0, 1]); bycov(ax, "fieldMode", "Field modal background (a.u.)")
plabel(ax, "B", x=-0.165)
ax = fig.add_subplot(gs[1, 0]); bycov(ax, "nCells", "Nuclei per field"); plabel(ax, "C")

ax = fig.add_subplot(gs[1, 1])
for f, g in order:
    v = np.sort(cell[cell.file == f].net.values)
    if v.size == 0:
        continue
    ax.step(v, np.arange(1, v.size + 1) / v.size, where="post", lw=0.9, color=GC[g],
            alpha=0.85, zorder=3)
ax.set_xscale("log"); ax.set_xlim(30, 4000); ax.set_ylim(-0.02, 1.03)
ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax.set_xlabel("Background-subtracted intensity per cell (a.u.)")
ax.set_ylabel("Cumulative fraction of cells")
ax.legend(handles=[Line2D([], [], color=GC[g], lw=1.2, label=FLAT[g]) for g in ALL],
          loc="upper left", frameon=False, handlelength=1.3, handletextpad=0.4,
          borderpad=0.05, labelspacing=0.22)
plabel(ax, "D", x=-0.165); tidy(ax)
save(fig, "FigureS12_pS6_fielddata")

# ---------------------------------------------------------------- numbers for the legends
print("\nS11 focus index by group (median):",
      {g: round(float(ok[ok.group == g].dapiNucEdge.median()), 3) for g in ALL})
print("S11 Kruskal-Wallis (4 treatment groups): P = %.4f" % kw)
print("S11 focus vs DAPI background: r = %+.3f, P = %.2e" % (R_FOCUS, P_FOCUS))
print("S11 P values by threshold:")
for c in CONTR:
    print("   %-18s %s" % ("%s vs %s" % (FLAT[c[0]], FLAT[c[1]]),
                           "  ".join("%.4f" % v for v in series[c])))
print("   thresholds", ths, " fields retained", kept)
print("\nS12 ring area median (um2):",
      {g: round(float(cell[cell.group == g].ringArea_um2.median()), 1) for g in ALL})
print("S12 pooled per-cell net vs ring area: r = %+.3f" %
      stats.pearsonr(cell.ringArea_um2, cell.net)[0])
print("\nS11a net vs nuclei per field (legend must match):")
for g in G:
    x = ok[ok.group == g]
    r, p = stats.pearsonr(x.nCells, x.net_field_true)
    print("   %-5s r = %+.2f, P = %.3f" % (FLAT[g], r, p))
