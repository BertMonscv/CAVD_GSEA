# Publication adaptation (2026-10-12): portable data/output path configuration only.
# Analysis, filtering, thresholds and plotting logic unchanged.
#!/usr/bin/env python3
# =============================================================================
#  Supplementary Figure S6 - VIC identity verification (vimentin / alpha-SMA)
#  CAVD manuscript - Molecular Medicine (BioMed Central) submission
#
#  FINAL production script. Run:  python make_FigureS6_MolMed.py
#  Outputs FigureS6_identity.{pdf,tif,png} plus a printed audit of every number
#  that appears in Supplementary Methods S2 and main-text Section 2.6, so the
#  figure and the manuscript text can be diffed against a single source.
#
#  ---------------------------------------------------------------------------
#  Journal specification applied (molmed.biomedcentral.com, Preparing figures)
#  ---------------------------------------------------------------------------
#   - full-page width exactly 170 mm; 107 mm tall, well inside the 225 mm limit
#   - closely cropped, minimal surrounding white space
#   - no figure title inside the graphic (titles and legends go in the text)
#   - figure KEYS are inside the graphic (threshold key in a/b, colour key in c)
#   - every rule >= 0.5 pt and every character >= 7 pt at final size
#   - fonts embedded as TrueType/Type 42, Arial/Helvetica metric family
#   - panel letters as bold lower-case a-d (BMC house style)
#   - TIFF with LZW compression and vector PDF, both far below the 10 MB limit
#
#  ---------------------------------------------------------------------------
#  Two deliberate design decisions, both defensible to a reviewer
#  ---------------------------------------------------------------------------
#  1. Central tendency is the MEDIAN with the INTERQUARTILE RANGE, matching
#     main-text Section 2.19 ("cell-level distributions ... summarized as
#     medians with interquartile ranges"). An earlier draft drew a median bar
#     with a mean +/- SD whisker, which mixed two summaries in one glyph.
#  2. Panels a and b have INDEPENDENT intensity axes. FITC exposure differed
#     between the two markers (100.2 ms for vimentin, 233.0 ms for alpha-SMA),
#     so raw counts are not comparable across markers and a shared axis would
#     imply a comparison the data do not support. Exposure was identical
#     between the stained and secondary-only coverslip WITHIN each marker -
#     the control that actually matters - and this is asserted at load time.
#     The legitimate cross-marker comparison is panel c, where each cell is
#     normalised to the mean of its own secondary-only control.
# =============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FixedFormatter, LogLocator, NullFormatter
from scipy.stats import gaussian_kde

from pathlib import Path
IDENTITY_DATA = Path(os.environ.get("IDENTITY_DATA", Path(__file__).resolve().parents[1] / "data"))
IDENTITY_FIGURES = Path(os.environ.get("IDENTITY_FIGURES", Path(__file__).resolve().parents[1] / "figures"))
IDENTITY_FIGURES.mkdir(parents=True, exist_ok=True)
PERFIELD_CSV = str(IDENTITY_DATA / "perfield.csv")
PERCELL_CSV = str(IDENTITY_DATA / "percell.csv")
OUTSTEM = str(IDENTITY_FIGURES / "FigureS6_identity")

# ---- house spec (shared with Figure 4D-G / Figure 5G) ----------------------
# Arial first; Liberation Sans is metric-compatible and used as a fallback so
# the layout is byte-identical on machines without Arial installed.
plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7.0,            # BMC-safe floor: nothing below 7 pt
    "axes.labelsize": 7.5,
    "xtick.labelsize": 7.0, "ytick.labelsize": 7.0,
    "legend.fontsize": 7.0,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
    "xtick.major.size": 2.6, "ytick.major.size": 2.6,
    "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
    "axes.grid": False, "legend.frameon": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
})

CTRL = "#AEB4B9"     # secondary-only control (house colour)
C_VIM = "#55606B"
C_SMA = "#8E44AD"
THR = "#A93226"      # threshold rules
REF = "#7F8C8D"      # 95 % reference rule

PT_S, PT_A = 1.4, 0.42       # per-cell point size / alpha
BEE_HW = 0.30                # beeswarm half-width
VIO_BW, VIO_HW, VIO_A = 0.35, 0.40, 0.13
BAR_HW, BAR_LW = 0.34, 1.35  # median bar
IQR_LW = 0.75                # interquartile whisker
SMA = "\u03b1-SMA"


# ---- primitives ------------------------------------------------------------
def bee(y, centre, hw=BEE_HW, nbin=None, log=True):
    """Deterministic binned jitter (no RNG), so the figure is reproducible."""
    v = np.log10(y) if log else np.asarray(y, float)
    if nbin is None:
        nbin = int(np.clip(4.5 * np.sqrt(len(v)), 40, 260))
    lo, hi = v.min(), v.max()
    if hi <= lo:
        hi = lo + 1e-9
    idx = np.clip(((v - lo) / (hi - lo) * (nbin - 1)).astype(int), 0, nbin - 1)
    x = np.empty(len(v))
    for b in np.unique(idx):
        s = np.where(idx == b)[0]
        s = s[np.argsort(v[s])]
        n = len(s)
        off = np.zeros(n) if n == 1 else np.linspace(-1, 1, n)
        x[s] = centre + off * hw * min(1.0, n / 7.0)
    return x


def violin(ax, y, centre, colour, log=True):
    if len(y) < 30:
        return
    v = np.log10(y) if log else np.asarray(y, float)
    k = gaussian_kde(v, bw_method=VIO_BW)
    gy = np.linspace(v.min(), v.max(), 220)
    d = k(gy)
    d = d / d.max() * VIO_HW
    ax.fill_betweenx(10 ** gy if log else gy, centre - d, centre + d,
                     color=colour, alpha=VIO_A, lw=0, zorder=1)


def tidy(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(which="both", direction="out")


def nice_log_ticks(lo, hi):
    """3-6 labelled ticks inside [lo, hi]: prefer 1-3-10, then 1-2-5, then decades."""
    dec = range(int(np.floor(np.log10(lo))), int(np.ceil(np.log10(hi))) + 1)
    best = None
    for mants in ([1, 3], [1, 2, 5], [1]):
        t = [mm * 10.0 ** e for e in dec for mm in mants]
        t = [x for x in t if lo <= x <= hi]
        if 3 <= len(t) <= 6:
            return t
        if best is None or abs(len(t) - 4) < abs(len(best) - 4):
            best = t
    return best


def plain_log_yaxis(ax, ticks):
    """Labelled ticks in plain numerals (BMC: no commas in numerical values)."""
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in ticks]))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=14))
    ax.yaxis.set_minor_formatter(NullFormatter())


def panel_letter(ax, letter, dx=-0.175):
    ax.annotate(letter, xy=(dx, 1.0), xycoords="axes fraction",
                xytext=(0, 9), textcoords="offset points",
                fontsize=8.5, fontweight="bold", va="top", ha="left",
                annotation_clip=False)


# ---- data ------------------------------------------------------------------
def load():
    """Field QC, then per-cell merge.

    Two inclusion criteria, applied identically to every coverslip:
      status == 'ok'   drops the computationally cleared duplicate that the
                       acquisition software writes for every field; only raw
                       16-bit frames are quantified.
      format match     keeps the standard 1024 x 1024 acquisition at 0.325 um
                       per pixel and drops fields acquired at 2048 x 2048.
                       This removes 2 of 24 secondary-only and 13 of 34
                       vimentin-stained fields, and nothing from alpha-SMA.
    """
    f = pd.read_csv(PERFIELD_CSV)
    c = pd.read_csv(PERCELL_CSV)
    ok = f[f.status == "ok"].copy()
    ok["matched"] = (ok.w == 1024) & (ok.px_um > 0.30)
    keep = ok[ok.matched]
    m = c.merge(keep[["file", "series", "w", "px_um", "expC1_s"]],
                on=["file", "series"], how="inner")

    D = {}
    for mk in ("VIM", "SMA"):
        d = m[m.marker == mk]
        ct, st = d[d.role == "control"], d[d.role == "stain"]
        exp_c = sorted(ct.expC1_s.unique())
        exp_s = sorted(st.expC1_s.unique())
        # The control that matters: identical FITC exposure within a marker.
        assert len(exp_c) == 1 and exp_c == exp_s, (
            f"{mk}: exposure differs between control {exp_c} and stain {exp_s}; "
            "intensities are not comparable and the thresholds are invalid.")
        raw = f[(f.marker == mk)]
        D[mk] = dict(
            ct=ct, st=st,
            T1=ct.ringMean.mean() + 3 * ct.ringMean.std(),
            T2=ct.ringMean.max(),
            cmean=ct.ringMean.mean(),
            exp_ms=exp_c[0] * 1000,
            fields_ok={r: int((raw.role == r).sum() // 2) for r in ("control", "stain")})
    return D


D = load()


def axis_range(mk):
    v = np.concatenate([D[mk]["ct"].ringMean.values, D[mk]["st"].ringMean.values])
    return v.min() / 1.28, v.max() * 1.28


# ---- figure ----------------------------------------------------------------
FIG_W_IN = 170 / 25.4        # BMC full page width
fig = plt.figure(figsize=(FIG_W_IN, 107 / 25.4))
outer = gridspec.GridSpec(2, 1, figure=fig, left=.077, right=.995, top=.950,
                          bottom=.130, hspace=.38)
top = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], wspace=.245)
bot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1], wspace=.245)


def percell_panel(ax, mk, colour, name):
    """Panels a and b: per-cell intensity, control vs stained coverslip."""
    d = D[mk]
    for i, (sub, col) in enumerate([(d["ct"], CTRL), (d["st"], colour)]):
        y = sub.ringMean.values
        violin(ax, y, i, col)
        ax.scatter(bee(y, i), y, s=PT_S, c=col, alpha=PT_A, linewidths=0,
                   rasterized=True, zorder=3)
        q1, med, q3 = np.percentile(y, [25, 50, 75])
        ax.vlines(i, q1, q3, color="k", lw=IQR_LW, zorder=4)
        ax.hlines(med, i - BAR_HW, i + BAR_HW, color="k", lw=BAR_LW, zorder=5)
    ax.axhline(d["T1"], color=THR, lw=.7, ls=(0, (4, 2)), zorder=2)
    ax.axhline(d["T2"], color=THR, lw=.7, ls=(0, (1, 1.6)), zorder=2)

    ax.set_yscale("log")
    lo, hi = axis_range(mk)
    ax.set_ylim(lo, hi)
    plain_log_yaxis(ax, nice_log_ticks(lo, hi))
    ax.set_xlim(-0.62, 1.62)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        [f"Secondary-only\n{len(d['ct'])} cells, {d['ct'].series.nunique()} fields",
         f"{name}\n{len(d['st'])} cells, {d['st'].series.nunique()} fields"])
    ax.set_ylabel("Perinuclear FITC intensity (a.u.)")
    tidy(ax)

    # threshold key inside the graphic (BMC: keys belong in the artwork)
    ax.legend(handles=[
        Line2D([], [], color=THR, lw=.7, ls=(0, (4, 2)),
               label=f"Control mean + 3 SD = {d['T1']:.0f}"),
        Line2D([], [], color=THR, lw=.7, ls=(0, (1, 1.6)),
               label=f"Control maximum = {d['T2']:.0f}")],
        loc="upper left", bbox_to_anchor=(0.015, 1.005),
        handlelength=2.1, handletextpad=0.5, labelspacing=0.30,
        borderpad=0.0, borderaxespad=0.0)


axA = fig.add_subplot(top[0])
percell_panel(axA, "VIM", C_VIM, "Vimentin")
axB = fig.add_subplot(top[1])
percell_panel(axB, "SMA", C_SMA, SMA)

# ---- c: fold over own control (normalises the exposure difference) ---------
axC = fig.add_subplot(bot[0])
grid = np.linspace(np.log10(0.25), np.log10(40), 400)
for mk, col, lab in [("VIM", C_VIM, "Vimentin"), ("SMA", C_SMA, SMA)]:
    d = D[mk]
    for sub, cc, ll, fill in [(d["ct"], CTRL, None, False), (d["st"], col, lab, True)]:
        v = np.log10(sub.ringMean.values / d["cmean"])
        k = gaussian_kde(v, bw_method=VIO_BW)(grid)
        if fill:
            axC.fill_between(10 ** grid, k, color=col, alpha=VIO_A, lw=0)
        axC.plot(10 ** grid, k, color=cc, lw=.9, label=ll)
axC.plot([], [], color=CTRL, lw=.9, label="Secondary-only controls")
axC.axvline(1, color="k", lw=.5, ls=(0, (4, 2)))
axC.set_xscale("log")
axC.set_xlim(0.25, 40)
axC.xaxis.set_major_locator(FixedLocator([0.3, 1, 3, 10, 30]))
axC.xaxis.set_major_formatter(FixedFormatter(["0.3", "1", "3", "10", "30"]))
axC.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=14))
axC.xaxis.set_minor_formatter(NullFormatter())
axC.tick_params(axis="x", which="minor", bottom=True)
axC.set_xlabel("Intensity / mean of own secondary-only control")
axC.set_ylabel("Density")
axC.legend(loc="upper right", handlelength=1.5, labelspacing=0.3, borderaxespad=0.2)
tidy(axC)

# ---- d: per-field positivity ----------------------------------------------
axD = fig.add_subplot(bot[1])
pos, ticks, labs, groups, audit = 0, [], [], [], {}
for mk, col, lab in [("VIM", C_VIM, "Vimentin"), ("SMA", C_SMA, SMA)]:
    d = D[mk]
    start = pos
    for T, tag in [(d["T1"], "mean + 3 SD"), (d["T2"], "control max")]:
        pf = d["st"].groupby("series").apply(
            lambda g: 100 * (g.ringMean > T).mean(), include_groups=False).values
        audit[(mk, tag)] = pf
        violin(axD, pf, pos, col, log=False)
        axD.scatter(bee(pf, pos, hw=0.24, nbin=30, log=False), pf, s=4.2, c=col,
                    alpha=.55, linewidths=0, rasterized=True, zorder=3)
        axD.hlines(pf.mean(), pos - BAR_HW, pos + BAR_HW, color="k", lw=BAR_LW, zorder=5)
        axD.text(pos, 100.9, f"{pf.mean():.1f}%", ha="center", va="bottom", fontsize=7)
        ticks.append(pos)
        labs.append(tag)
        pos += 1
    groups.append((lab, start, pos - 1))

axD.axhline(95, color=REF, lw=.6, ls=(0, (4, 2)), zorder=2)
axD.text(pos - 0.45, 95, "95%", color=REF, fontsize=7, ha="right", va="bottom")
axD.set_xticks(ticks)
axD.set_xticklabels(labs)
axD.set_ylim(84, 102.4)
axD.set_yticks([85, 90, 95, 100])
axD.set_xlim(-0.62, pos - 0.38)
axD.set_ylabel("Positive cells per field (%)")
axD.spines["left"].set_bounds(85, 100)   # axis truncated at 85 %, stated in the legend
tidy(axD)
axD.minorticks_on()
axD.tick_params(axis="x", which="minor", bottom=False)
axD.tick_params(axis="y", which="minor", left=False)

# marker group brackets under the tick labels
tr = axD.get_xaxis_transform()
for lab, a, b in groups:
    axD.plot([a - 0.42, b + 0.42], [-0.175, -0.175], transform=tr, color="k",
             lw=0.6, clip_on=False)
    axD.text((a + b) / 2, -0.215, lab, transform=tr, ha="center", va="top",
             fontsize=7.5, clip_on=False)

for ax, L in [(axA, "a"), (axB, "b"), (axC, "c"), (axD, "d")]:
    panel_letter(ax, L)

# ---- output ----------------------------------------------------------------
# No bbox_inches: the canvas is already exactly 170 mm and closely cropped;
# tight cropping would silently shrink the width to ~164 mm.
fig.savefig(OUTSTEM + ".pdf")                                             # vector
fig.savefig(OUTSTEM + ".tif", dpi=400, pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(OUTSTEM + ".png", dpi=600)                                    # screen check

# Matplotlib writes RGBA; flatten the TIFF onto white and drop the alpha channel,
# because several production pipelines reject TIFFs carrying transparency.
from PIL import Image

with Image.open(OUTSTEM + ".tif") as _im:
    _rgba = _im.convert("RGBA")
    _flat = Image.alpha_composite(Image.new("RGBA", _rgba.size, (255, 255, 255, 255)),
                                  _rgba).convert("RGB")
    _flat.save(OUTSTEM + ".tif", compression="tiff_lzw", dpi=(400, 400))

# ---- compliance and numbers audit ------------------------------------------
w_in, h_in = fig.get_size_inches()
print("=" * 74)
print(f"width  {w_in * 25.4:6.1f} mm   (BMC full page 170 mm)")
print(f"height {h_in * 25.4:6.1f} mm   (BMC ceiling  225 mm)")
for ext in (".pdf", ".tif", ".png"):
    p = OUTSTEM + ext
    print(f"{ext:5s}  {os.path.getsize(p) / 1e6:5.2f} MB   (BMC ceiling 10 MB)")
print("=" * 74)
print("Numbers below must match Supplementary Methods S2 and Section 2.6.\n")
for mk, name in [("VIM", "Vimentin"), ("SMA", SMA)]:
    d = D[mk]
    ct, st = d["ct"].ringMean.values, d["st"].ringMean.values
    print(f"--- {name}  (FITC exposure {d['exp_ms']:.1f} ms, identical for both coverslips)")
    print(f"    fields acquired -> retained : control {d['fields_ok']['control']} -> "
          f"{d['ct'].series.nunique()} | stain {d['fields_ok']['stain']} -> {d['st'].series.nunique()}")
    print(f"    control  n = {len(ct):4d}   {ct.mean():7.1f} +/- {ct.std(ddof=1):5.1f}   "
          f"p99 = {np.percentile(ct, 99):.0f}   max = {ct.max():.0f}")
    print(f"    stained  n = {len(st):4d}   {st.mean():7.1f} +/- {st.std(ddof=1):5.1f}   "
          f"median = {np.median(st):.0f} ({np.median(st) / d['cmean']:.1f}x control mean)   "
          f"IQR {np.percentile(st, 25):.0f}-{np.percentile(st, 75):.0f}   min = {st.min():.0f}")
    for T, tag in [(d["T1"], "mean + 3 SD"), (d["T2"], "control max")]:
        pf = audit[(mk, tag)]
        print(f"    {tag:12s} = {T:4.0f} : stained {int((st > T).sum())}/{len(st)} "
              f"({100 * (st > T).mean():6.2f}%) | control false-positive "
              f"{int((ct > T).sum())}/{len(ct)} ({100 * (ct > T).mean():.2f}%) | "
              f"per-field {pf.mean():.2f} +/- {pf.std(ddof=1):.2f}%, lowest {pf.min():.1f}%")
    print()
