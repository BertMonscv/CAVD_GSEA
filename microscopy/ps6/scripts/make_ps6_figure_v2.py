# Publication copy (2026-10-12): analysis and plotting logic unchanged.
# Historical output figure labels are mapped to current manuscript figures in README.
#!/usr/bin/env python3
"""Figure 4C-F — p-S6 (Ser235/236) single-cell quantification. Journal-submission build.

Changes over the first pass:
  * C: violin envelope behind the beeswarm so all four groups read at a glance;
       mean +/- SD whisker on the coverslip medians; n printed under each group.
  * D: ridgeline KDEs replace the jagged step histograms (the 2-only group had
       only 76 cells and looked like noise); median tick per group.
  * E: curve labels at the line ends instead of a legend box and a crossing arrow.
  * F: per-group values printed; paired-t comparison against CTL added.
  * Global: colour-blind-safe palette, outward ticks, minor ticks, hairline spines,
       no chartjunk, per-cell points rasterised, fonts embedded as Type 42.
Run with no arguments for the current Figure 4 (panels a-d). PANEL_SET can be set to
"DEFG" or "CDEF" to regenerate the superseded panel-lettering used before Figure 3 was split.
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator
from scipy import stats
from scipy.stats import gaussian_kde

PANEL_SET = os.environ.get("PANEL_SET", "ABCD")
LABELS = list(PANEL_SET)
STEMS = {"ABCD": "Figure4_pS6_quant",       # current: p-S6 quantification is Figure 4
         "DEFG": "Figure4DG_pS6_quant",   # superseded, kept for reference
         "CDEF": "Figure4CF_pS6_quant"}   # superseded
OUTSTEM = STEMS.get(PANEL_SET, "Figure_pS6_quant")

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
cell["rep"] = cell["rep"].astype(str)

G = ["CTL", "PM", "LEP", "RAPA"]
LAB = {"CTL": "CTL", "PM": "PM", "LEP": "PM+Lep", "RAPA": "PM+Lep\n+Rapa", "Sac": "2°-only"}
FLAT = {k: v.replace("\n", " ") for k, v in LAB.items()}
GC = {"CTL": "#55606B", "PM": "#2E6DB4", "LEP": "#C0392B", "RAPA": "#1E8A5A", "Sac": "#AEB4B9"}
RC = {"1": "#4E79A7", "2": "#F28E2B", "3": "#59A14F"}

fieldmed = (cell[cell.group.isin(G)]
            .groupby(["group", "rep", "file", "series"], as_index=False).net.median())
piv = (fieldmed.groupby(["group", "rep"], as_index=False).net.median()
       .pivot(index="rep", columns="group", values="net")[G])

CONTR = [("CTL", "PM"), ("PM", "LEP"), ("LEP", "RAPA")]
praw = {c: stats.ttest_rel(piv[c[0]], piv[c[1]]).pvalue for c in CONTR}
order = sorted(praw, key=praw.get)
padj = {c: min(1.0, praw[c] * (len(order) - i)) for i, c in enumerate(order)}
for i in range(1, len(order)):
    padj[order[i]] = max(padj[order[i]], padj[order[i - 1]])

sac = cell[cell.group == "Sac"].net
ctl95, sac95 = cell[cell.group == "CTL"].net.quantile(.95), sac.quantile(.95)
NCELL = {g: int((cell.group == g).sum()) for g in G + ["Sac"]}
NFIELD = fieldmed.groupby("group").size().to_dict()


def pf(p):
    return "P < 0.001" if p < 0.001 else "P = %.3f" % p


def swarm(y, cx, halfwidth, rng, nbin=40):
    y = np.asarray(y, float)
    x = np.full(y.size, cx)
    lo, hi = y.min(), y.max()
    if hi <= lo:
        return x
    b = np.clip(((y - lo) / (hi - lo) * (nbin - 1)).astype(int), 0, nbin - 1)
    for k in range(nbin):
        m = np.where(b == k)[0]
        if not m.size:
            continue
        off = np.arange(m.size) - (m.size - 1) / 2
        span = min(halfwidth, halfwidth * m.size / 13.0)
        step = 0 if m.size < 2 else span / (m.size - 1)
        x[m] = cx + off * step + rng.uniform(-.005, .005, m.size)
    return x


def plabel(ax, ch, x=-0.185, y=1.048):
    ax.text(x, y, ch.lower(), transform=ax.transAxes, fontsize=9.5,
            fontweight="bold", ha="left", va="bottom")


def tidy(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(which="both", pad=1.8)


rng = np.random.default_rng(20260804)
fig = plt.figure(figsize=(6.69, 5.90))
gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.30,
                      left=0.112, right=0.984, top=0.945, bottom=0.090)

# ══════════════════════════════════ panel 1 — SuperPlot ══════════════════════
ax = fig.add_subplot(gs[0, 0])
ax.axhspan(np.percentile(sac, 25), np.percentile(sac, 75), color="#E4E6E8", lw=0, zorder=0)
ax.axhline(sac.median(), color="#9AA0A6", lw=.65, ls=(0, (4, 2.5)), zorder=1)

for i, g in enumerate(G):
    v = cell[cell.group == g].net.values
    kde = gaussian_kde(v, bw_method=.35)
    yy = np.linspace(v.min(), v.max(), 220)
    w = kde(yy); w = w / w.max() * .40
    ax.fill_betweenx(yy, i - w, i + w, facecolor=GC[g], alpha=.13, lw=0, zorder=1)
    ax.plot(i - w, yy, color=GC[g], lw=.55, zorder=2)
    ax.plot(i + w, yy, color=GC[g], lw=.55, zorder=2)
    for r in "123":
        s = cell[(cell.group == g) & (cell.rep == r)].net.values
        ax.scatter(swarm(s, i + (int(r) - 2) * .205, .145, rng), s,
                   s=1.7, c=RC[r], alpha=.42, lw=0, rasterized=True, zorder=3)
    for r in "123":
        ax.scatter(i + (int(r) - 2) * .205, piv.loc[r, g], s=34, marker="D",
                   facecolor=RC[r], edgecolor="black", lw=.6, zorder=6)
    m, sd = piv[g].mean(), piv[g].std(ddof=1)
    ax.plot([i - .36, i + .36], [m, m], color="black", lw=1.5, zorder=7,
            solid_capstyle="butt")
    ax.plot([i, i], [m - sd, m + sd], color="black", lw=.75, zorder=5)

for (a, b), h in zip(CONTR, [1380, 3330, 3610]):
    ia, ib = G.index(a), G.index(b)
    ax.plot([ia, ia, ib, ib], [h, h + 65, h + 65, h], lw=.55, c="black", clip_on=False)
    ax.text((ia + ib) / 2, h + 90, pf(padj[(a, b)]), ha="center", va="bottom", fontsize=7.4)

ax.set_xticks(range(4)); ax.set_xticklabels([LAB[g] for g in G])
ax.set_xlim(-.62, 3.62); ax.set_ylim(-70, 3980)
for i, g in enumerate(G):
    ax.text(i, -0.20, "%d / %d" % (NCELL[g], NFIELD[g]), transform=
            ax.get_xaxis_transform(), ha="center", va="top", fontsize=6.9, color="#4a4a4a")
ax.annotate("2°-only control\n(median, IQR)", xy=(3.52, sac.median()), xytext=(3.52, 1150),
            fontsize=6, ha="right", va="bottom", color="#6b7075",
            arrowprops=dict(arrowstyle="-", lw=.45, color="#9AA0A6"))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.set_ylabel("Background-subtracted p-S6\nintensity per cell (a.u.)")
plabel(ax, LABELS[0])
tidy(ax)

# ══════════════════════════════════ panel 2 — ridgeline ══════════════════════
ax = fig.add_subplot(gs[0, 1])
rows = ["LEP", "PM", "CTL", "RAPA", "Sac"]
grid = np.linspace(1.35, 3.62, 400)
for j, g in enumerate(rows):
    v = np.log10(cell[cell.group == g].net.clip(lower=1))
    d = gaussian_kde(v, bw_method=.30)(grid)
    d = d / d.max() * .84
    ax.fill_between(grid, j, j + d, facecolor=GC[g], alpha=.55, lw=0, zorder=3 + j)
    ax.plot(grid, j + d, color=GC[g], lw=.9, zorder=3 + j)
    ax.plot([grid[0], grid[-1]], [j, j], color="#c9cdd1", lw=.4, zorder=2)
    med = np.median(v)
    ax.plot([med, med], [j, j + np.interp(med, grid, d)], color="white", lw=.9, zorder=9 + j)
    ax.text(1.33, j + .12, FLAT[g], ha="right", va="bottom", fontsize=7.6, color=GC[g])
ax.axvline(np.log10(ctl95), color="black", lw=.6, ls=":", zorder=20)
ax.text(np.log10(ctl95) + .03, 4.95, "CTL 95th\npercentile", fontsize=6.9, va="top")
ax.set_xlim(1.05, 3.62); ax.set_ylim(-.15, 5.05)
ax.set_yticks([])
ax.set_xticks([1.5, 2.0, 2.5, 3.0, 3.5])
ax.set_xlabel("p-S6 intensity per cell (log$_{10}$ a.u.)")
plabel(ax, LABELS[1], x=-0.045)
ax.spines[["top", "right", "left"]].set_visible(False)
ax.tick_params(which="both", pad=1.8)

# ══════════════════════════════════ panel 3 — ECDF ═══════════════════════════
ax = fig.add_subplot(gs[1, 0])
ends = {"Sac": .875, "RAPA": .815, "CTL": .755, "PM": .695, "LEP": .635}
for g in ["Sac", "RAPA", "CTL", "PM", "LEP"]:
    v = np.sort(cell[cell.group == g].net.values)
    ax.step(v, np.arange(1, v.size + 1) / v.size, where="post", lw=1.25, color=GC[g], zorder=4)
    ax.text(3350, ends[g], FLAT[g], color=GC[g], fontsize=7.4, ha="right", va="center")
ax.axvline(ctl95, color="black", lw=.6, ls=":", zorder=2)
ax.axvline(sac95, color="#9AA0A6", lw=.6, ls=(0, (4, 2.5)), zorder=2)
ax.set_xscale("log"); ax.set_xlim(33, 5200); ax.set_ylim(-.02, 1.03)
ax.set_yticks([0, .25, .5, .75, 1])
frac = 100 * (cell[cell.group == "LEP"].net > ctl95).mean()
ax.plot([ctl95], [1 - frac / 100], marker="o", ms=3.2, mfc="white",
        mec=GC["LEP"], mew=.9, zorder=6)
ax.annotate("%.1f%% of PM+Lep cells\nabove CTL 95th pct" % frac,
            xy=(ctl95 * 1.07, 1 - frac / 100), xytext=(5000, .045), fontsize=7.2,
            va="bottom", ha="right",
            arrowprops=dict(arrowstyle="-", lw=.45, color="#8a8f94", shrinkB=3))
ax.text(ctl95, 1.035, "CTL 95th percentile", fontsize=6.9, ha="center", va="bottom",
        color="#4a4a4a", clip_on=False)
ax.set_xlabel("p-S6 intensity per cell (a.u., log scale)")
ax.set_ylabel("Cumulative fraction of cells")
plabel(ax, LABELS[2])
tidy(ax)

# ══════════════════════════════════ panel 4 — fraction above ═════════════════
ax = fig.add_subplot(gs[1, 1])
fr = {g: [100 * (cell[(cell.group == g) & (cell.rep == r)].net > ctl95).mean() for r in "123"]
      for g in G}
for i, g in enumerate(G):
    for r in "123":
        ax.scatter(i + (int(r) - 2) * .175, fr[g][int(r) - 1], s=31, marker="D",
                   facecolor=RC[r], edgecolor="black", lw=.5, zorder=5)
    m, sd = np.mean(fr[g]), np.std(fr[g], ddof=1)
    ax.plot([i - .33, i + .33], [m, m], color="black", lw=1.45, zorder=6, solid_capstyle="butt")
    ax.plot([i, i], [m - sd, m + sd], color="black", lw=.75, zorder=4)
    ax.text(i, min(m + sd + 4.5, 96), "%.1f" % m, ha="center", va="bottom", fontsize=7.3)
ax.axhline(5, color="#9AA0A6", lw=.6, ls=":", zorder=1)
ax.text(3.58, 15, "5% expected\nby definition", fontsize=6.8, ha="right", va="bottom", color="#6b7075")
ax.set_xticks(range(4)); ax.set_xticklabels([LAB[g] for g in G])
ax.set_xlim(-.62, 3.62); ax.set_ylim(-5, 103)
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.set_ylabel("Cells above CTL 95th percentile (%)")
plabel(ax, LABELS[3])
tidy(ax)
ax.legend(handles=[Line2D([], [], marker="D", ls="", mfc=RC[r], mec="k", mew=.5, ms=4.4,
                          label="Experiment %s" % r) for r in "123"],
          loc="center left", bbox_to_anchor=(0.0, 0.62), frameon=False,
          handletextpad=.35, borderpad=.05, labelspacing=.28)

for ext in ("png", "pdf"):
    fig.savefig(OUTDIR / ("%s.%s" % (OUTSTEM, ext)), dpi=600,
                facecolor="white")   # no bbox trim: canvas is exactly 170 mm wide
plt.close(fig)
from PIL import Image
Image.open(OUTDIR / ("%s.png" % OUTSTEM)).convert("RGB").save(
    OUTDIR / ("%s.tif" % OUTSTEM), compression="tiff_lzw", dpi=(600, 600))
print("wrote", OUTSTEM, "panels", LABELS)
