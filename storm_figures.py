"""STORM synapse figures — Python implementation mirroring the MATLAB code.

Provides:
    load_storm_data(path)        -> dict of pandas DataFrames
    plot_density(data, age)      -> matplotlib Figure  (Fig 4B / S6A style)
    plot_volume_intensity(...)   -> matplotlib Figure  (Fig 4D / S6C style)
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

COLORS = {"Het": (0, 0, 0), "KO": (0.85, 0.10, 0.10)}
# Beeswarm point-cloud colours: Het/WT in gray, KO in red.
CLOUD_COLORS = {"Het": (0.55, 0.55, 0.55), "KO": (0.85, 0.10, 0.10)}
BLACK = (0, 0, 0)
GENOS = ["Het", "KO"]


def _dodge_overlaps(y, min_gap, step):
    """Return x-offsets that spread points colliding in y around the center.

    Points whose y-values fall within ``min_gap`` of a neighbour are treated
    as a cluster and fanned out symmetrically about x = 0 with spacing
    ``step`` (…, -step, 0, +step, … for odd counts; ±step/2, ±3step/2, … for
    even counts). Isolated points stay centred at 0.
    """
    y = np.asarray(y, dtype=float)
    n = y.size
    offs = np.zeros(n)
    if n == 0:
        return offs
    order = np.argsort(y)
    ys = y[order]
    # group consecutive points that are within min_gap of the previous one
    group_id = np.zeros(n, dtype=int)
    gid = 0
    for i in range(1, n):
        if ys[i] - ys[i - 1] < min_gap:
            group_id[i] = gid
        else:
            gid += 1
            group_id[i] = gid
    for g in np.unique(group_id):
        idx = np.where(group_id == g)[0]
        m = idx.size
        if m == 1:
            continue
        # symmetric positions centred on 0
        pos = (np.arange(m) - (m - 1) / 2.0) * step
        offs[order[idx]] = pos
    return offs

_VI_SHEETS = {
    "VGluT2": "VGluT2_volume_intensity",
    "Bassoon": "Bassoon_volume_intensity",
    "Homer": "Homer_volume_intensity",
}


def _norm_age(s):
    s = str(s).strip().lower()
    if s in ("p8", "p08"):
        return "P8"
    if s in ("adult", "p60", "p 60"):
        return "Adult"
    return str(s)


def _norm_geno(s):
    s = str(s).strip().lower()
    if s in ("het", "opn4-het", "opn4cre/+"):
        return "Het"
    if s in ("ko", "opn4-ko", "opn4cre/cre"):
        return "KO"
    return str(s)


def load_storm_data(path=None):
    """Read the workbook into a dict of tidy DataFrames."""
    if path is None:
        path = Path(__file__).resolve().parent / "data" / "STORM_raw_spreadsheet.xlsx"
    path = Path(path)
    xl = pd.ExcelFile(path)

    dens = pd.read_excel(xl, "Synapse_density")
    dens = dens[["Sample_ID", "Age", "Genotype",
                 "Biological replicate", "Channel", "Density"]].copy()
    dens.columns = ["Sample_ID", "Age", "Genotype", "Replicate", "Channel", "Density"]
    dens["Age"] = dens["Age"].map(_norm_age)
    dens["Genotype"] = dens["Genotype"].map(_norm_geno)
    dens = dens.dropna(subset=["Density"])

    data = {"density": dens}
    for prot, sheet in _VI_SHEETS.items():
        df = pd.read_excel(xl, sheet)
        df = df[["Sample_ID", "Age", "Genotype", "Biological replicate",
                 "Volume log(um3)", "Intensity log(A.U.)"]].copy()
        df.columns = ["Sample_ID", "Age", "Genotype", "Replicate",
                      "Volume", "Intensity"]
        df["Age"] = df["Age"].map(_norm_age)
        df["Genotype"] = df["Genotype"].map(_norm_geno)
        df = df.dropna(subset=["Volume", "Intensity"])
        data[prot] = df
    return data


def plot_density(data, age, channel="VGluT2", ax=None):
    """Dot plot of per-replicate synapse density with mean +/- S.D.

    Individual animal values are aligned centrally along the error bar
    (no horizontal jitter)."""
    df = data["density"]
    df = df[(df["Age"] == age) & (df["Channel"] == channel)]

    if ax is None:
        fig, ax = plt.subplots(figsize=(2.2, 3.0))
    else:
        fig = ax.figure

    xpos = {"Het": 1, "KO": 2}
    yrange = df["Density"].max() * 1.15  # matches ylim below
    for g in GENOS:
        v = df.loc[df["Genotype"] == g, "Density"].to_numpy()
        c = COLORS[g]
        mu, sd = v.mean(), v.std(ddof=1)
        # individual animals near the error bar; nudge apart if they overlap
        dx = _dodge_overlaps(v, min_gap=0.03 * yrange, step=0.055)
        ax.scatter(xpos[g] + dx, v, s=30, color=c, alpha=0.9,
                   edgecolors=c, zorder=3)
        ax.errorbar(xpos[g], mu, yerr=sd, color=c, lw=1.2, capsize=8,
                    capthick=1.2, zorder=4)
        ax.plot([xpos[g] - 0.22, xpos[g] + 0.22], [mu, mu],
                color=c, lw=1.6, zorder=5)

    ax.set_xlim(0.4, 2.6)
    ax.set_xticks([1, 2]); ax.set_xticklabels(["Het", "KO"])
    ax.set_ylabel(r"Synapse density (per $\mu$m$^3$)")
    ax.set_ylim(0, df["Density"].max() * 1.15)
    ax.set_title(f"{channel}  ({age})")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out")
    fig.tight_layout()
    return fig


def _violin_offsets(v, width=0.36, rng_seed=0):
    """Horizontal offset for each point so the cloud forms a violin.

    Every synapse is drawn as a circle. The half-width available to a point
    at value y is proportional to the local density there (a KDE evaluated at
    the point), and each point is placed at a random position within that
    half-width. The envelope of the resulting cloud is the violin silhouette,
    which is generated by the points themselves rather than drawn as a curve.
    """
    v = np.asarray(v, dtype=float)
    n = v.size
    if n == 0:
        return np.array([]), v
    if n < 2 or np.allclose(v, v[0]):
        return np.zeros(n), v
    dens = gaussian_kde(v)(v)
    dens = dens / dens.max()               # 0..1 local density
    rng = np.random.default_rng(rng_seed)
    off = (rng.random(n) * 2 - 1) * dens * width
    return off, v


def _violin_box(ax, xc, v, cloud_c, animal_means=None, rng_seed=0, yspan=None):
    v = np.asarray(v)
    v = v[~np.isnan(v)]
    if v.size < 2:
        return

    # --- point cloud: one circle per synapse, shaped as a violin ---
    # (Het/WT in gray, KO in red)
    off, yy = _violin_offsets(v, width=0.36, rng_seed=rng_seed)
    ax.scatter(xc + off, yy, s=1.5, color=cloud_c, alpha=0.25,
               edgecolors="none", rasterized=True, zorder=1)

    # --- per-animal means (N=6) overlaid centrally as black circles ---
    if animal_means is not None and len(animal_means) > 0:
        am = np.asarray(animal_means, dtype=float)
        am = am[~np.isnan(am)]
        yr = yspan if yspan else (am.max() - am.min() + 1e-9)
        dx = _dodge_overlaps(am, min_gap=0.035 * yr, step=0.05)
        ax.scatter(xc + dx, am, s=16, color=BLACK, alpha=0.95,
                   edgecolors="white", linewidths=0.4, zorder=6)

    # --- box-and-whisker overlay (always black lines) ---
    q1, q2, q3 = np.percentile(v, [25, 50, 75])
    iqr = q3 - q1
    lo = max(v.min(), q1 - 1.5 * iqr)
    hi = min(v.max(), q3 + 1.5 * iqr)
    bw = 0.12
    ax.plot([xc, xc], [lo, q1], color=BLACK, lw=1.0, zorder=3)
    ax.plot([xc, xc], [q3, hi], color=BLACK, lw=1.0, zorder=3)
    ax.plot([xc - bw, xc + bw], [lo, lo], color=BLACK, lw=1.0, zorder=3)
    ax.plot([xc - bw, xc + bw], [hi, hi], color=BLACK, lw=1.0, zorder=3)
    ax.add_patch(plt.Rectangle((xc - bw, q1), 2 * bw, q3 - q1,
                               edgecolor=BLACK, facecolor="none", lw=1.2, zorder=4))
    ax.plot([xc - bw, xc + bw], [q2, q2], color=BLACK, lw=1.6, zorder=5)


def plot_volume_intensity(data, age):
    """Six-panel violin+box figure: Volume & Intensity x 3 proteins."""
    proteins = ["VGluT2", "Bassoon", "Homer"]
    prot_lbl = {"VGluT2": "VGluT2", "Bassoon": "Bassoon", "Homer": "Homer1"}
    measures = ["Volume", "Intensity"]
    meas_lbl = {"Volume": r"Volume log($\mu$m$^3$)",
                "Intensity": "Intensity log(A.U.)"}

    fig, axes = plt.subplots(1, 6, figsize=(11, 3.4))
    col = 0
    for prot in proteins:
        df = data[prot]
        df = df[df["Age"] == age]
        for meas in measures:
            ax = axes[col]; col += 1
            panel_span = df[meas].max() - df[meas].min()
            for i, g in enumerate(GENOS, start=1):
                sub = df[df["Genotype"] == g]
                v = sub[meas].to_numpy()
                means = sub.groupby("Replicate")[meas].mean().to_numpy()
                _violin_box(ax, i, v, CLOUD_COLORS[g], animal_means=means,
                            yspan=panel_span)
            ax.set_xlim(0.4, 2.6)
            ax.set_xticks([1, 2]); ax.set_xticklabels(["Het", "KO"])
            ax.set_ylabel(meas_lbl[meas])
            ax.set_title(prot_lbl[prot])
            ax.spines[["top", "right"]].set_visible(False)
            ax.tick_params(direction="out")
    fig.suptitle(f"{age} synaptic protein distributions", y=1.02)
    fig.tight_layout()
    return fig
