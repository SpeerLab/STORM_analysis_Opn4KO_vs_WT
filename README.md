# STORM RHT Synapse Figures

Code to reproduce the STORM synapse-analysis figures for the melanopsin / RHT
synaptogenesis study:

- **Figure 4** — P8 (prior to eye-opening)
- **Figure S6** — Adult

Two parallel implementations are provided so results can be regenerated in
either language, and a Jupyter notebook renders everything inline for viewing
directly on GitHub:

- `matlab/` — MATLAB functions + a one-shot driver script
- `storm_figures.py` — the same plots in Python (matplotlib + SciPy)
- `STORM_figures.ipynb` — executed notebook with the figures embedded

## What gets plotted

**Synapse density** (`plot_density`) — one dot per biological replicate
(N = 6 per genotype), with the group **mean ± S.D.** overlaid as an error bar.
Het = black, KO = red.

**Volume & intensity** (`plot_volume_intensity`) — violin + box-and-whisker
panels for **VGluT2**, **Bassoon**, and **Homer1**, showing volume (log µm³)
and intensity (log A.U.). The box marks the median and interquartile range,
whiskers extend to 1.5× IQR, and the violin outline is a kernel-density
estimate of the full distribution.

## Repository layout

```
STORM_figures/
├── data/
│   └── STORM_raw_spreadsheet.xlsx     # source data (4 sheets)
├── matlab/
│   ├── load_storm_data.m              # workbook -> struct of tables
│   ├── plot_density.m                 # density dot plot
│   ├── plot_volume_intensity.m        # 6-panel violin + box figure
│   └── make_all_figures.m             # runs everything, writes to output/
├── storm_figures.py                   # Python equivalent of the above
├── STORM_figures.ipynb                # executed notebook (view on GitHub)
├── output/                            # generated .eps / .png figures
├── requirements.txt
└── README.md
```

## Data

`data/STORM_raw_spreadsheet.xlsx` contains four sheets:

| Sheet | Rows | Contents |
|---|---|---|
| `Synapse_density` | 1 per replicate | density per genotype/age/channel |
| `VGluT2_volume_intensity` | 1 per cluster | volume & intensity, VGluT2 |
| `Bassoon_volume_intensity` | 1 per cluster | volume & intensity, Bassoon |
| `Homer_volume_intensity` | 1 per cluster | volume & intensity, Homer1 |

Both ages (`P8`, `Adult`) and both genotypes (`Het`, `KO`) live in every sheet;
the plotting functions filter by age.

## Running — Python

```bash
pip install -r requirements.txt
jupyter notebook STORM_figures.ipynb        # interactive
# or headless:
jupyter nbconvert --to notebook --execute --inplace STORM_figures.ipynb
```

Or from a script / REPL:

```python
import storm_figures as sf
data = sf.load_storm_data("data/STORM_raw_spreadsheet.xlsx")
sf.plot_density(data, "Adult")            # Figure S6A
sf.plot_volume_intensity(data, "Adult")   # Figure S6C
sf.plot_density(data, "P8")               # Figure 4B
sf.plot_volume_intensity(data, "P8")      # Figure 4D
```

## Running — MATLAB

From the `matlab/` directory (or with it on the path):

```matlab
make_all_figures        % writes .eps + .png for both ages into ../output/
```

Or call the functions directly:

```matlab
data = load_storm_data();                                  % ../data/…xlsx
plot_density(data, 'Adult', 'SaveStem', '../output/S6A');  % Figure S6A
plot_volume_intensity(data, 'Adult', 'SaveStem', '../output/S6C');
```

MATLAB R2020a+ is recommended (`tiledlayout`, `exportgraphics`).
`ksdensity` requires the Statistics and Machine Learning Toolbox.

## Outputs

Each figure is saved as a vector **`.pdf`** plus a 300-dpi **`.png`**:

- `Fig4B_density_P8`, `Fig4D_volint_P8`
- `FigS6A_density_Adult`, `FigS6C_volint_Adult`

**Why PDF, not EPS?** The violin panels draw every synapse as a
semi-transparent circle, and that transparency is what gives the point cloud
its soft, density-graded shape. The EPS format cannot store transparency —
exporting to EPS flattens the points to fully opaque and the clouds collapse
into solid blobs. PDF keeps the transparency and is fully vector, and
Illustrator opens it natively. If a journal specifically requires `.eps`,
open the `.pdf` in Illustrator and *Save As → EPS*; Illustrator rasterises the
transparent layer correctly on export, which the plotting libraries cannot.

## How the plots are built

**Density** — the six per-animal values are drawn as dots aligned centrally on
the error bar (no horizontal jitter), with the mean line and ±1 S.D. whiskers.

**Volume / intensity violins** — the violin shape is **not** a drawn kernel
outline. Every synapse is plotted as its own circle, and each point's
horizontal offset is scaled by the local point density (so crowded values
spread wider). The envelope of the point cloud *is* the violin. The larger
opaque dots overlaid centrally are the six per-animal means; the
box-and-whisker (median, quartiles, 1.5× IQR whiskers) sits on top.

## Notes

- The mean ± S.D. values reproduce those quoted in the figure legends
  (e.g. Adult Het 0.044 ± 0.007, KO 0.045 ± 0.003; P8 Het 0.030 ± 0.004,
  KO 0.023 ± 0.002).
- Statistical testing is **not** performed by this code; the plots visualise
  the distributions only. The p-values / Cohen's d in the legends should be
  computed and reported separately.
- The point-cloud offsets use a fixed random seed so the violins are
  reproducible run to run.
