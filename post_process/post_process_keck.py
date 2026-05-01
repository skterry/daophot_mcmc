import sys
import os
import getdist
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image, ImageOps, ImageChops
from astropy.visualization import astropy_mpl_style
import corner

from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

import pdb

import smplotlib
#plt.style.use(astropy_mpl_style)

#plt.rcParams.update({'font.size': 20})
plt.rcParams["font.family"] = "courier new"

# Parse star-count argument: e.g. "python post_process_keck.py -2"
if len(sys.argv) != 2 or not (sys.argv[1].startswith('-') and sys.argv[1][1:].isdigit()):
    print("Usage: python post_process_keck.py -<n>  (e.g., -2 for a 2-star fit)")
    sys.exit(1)

n_arg = int(sys.argv[1][1:])
fname   = f'../{n_arg}star.mcmc'
c_fname = f'../{n_arg}star.chi2pix'
bf_fname = f'../{n_arg}star.bestfit'

for fpath in [fname, c_fname]:
    if not os.path.exists(fpath):
        print(f"Cannot run: file '{fpath}' does not exist.")
        sys.exit(1)

#New format read FITS and display pixels/values:
image_data = fits.getdata('../image.fits', ext=0)

with open(fname) as f:
    col_names = f.readline().lstrip('#').split()
chains_all = pd.read_table(fname, sep=r'\s+', skiprows=1, names=col_names)

# Determine number of stars from column count:
#   1-star: X, Y, Flux, chi^2                           (4 cols)
#   2-star: X1, Y1, X2, Y2, Flux1, Flux2, Sep, chi^2   (8 cols)
#   3-star: X1, Y1, X2, Y2, X3, Y3, ...                (>8 cols)
n_cols = len(col_names)
if n_cols <= 4:
    n_stars = 1
elif n_cols <= 8:
    n_stars = 2
else:
    n_stars = 3

star_suffix = f'{n_stars}star'

# Extract per-star XY chains with Python index correction (-1)
star_xy = []
for i in range(n_stars):
    xcol, ycol = col_names[2 * i], col_names[2 * i + 1]
    sc = chains_all[[xcol, ycol]].copy()
    sc[xcol] -= 1
    sc[ycol] -= 1
    star_xy.append((xcol, ycol, sc))

# Remaining non-XY columns (flux, separation, chi^2, etc.)
extra_col_names = col_names[2 * n_stars:]
chains_extra = chains_all[extra_col_names].copy() if extra_col_names else pd.DataFrame()

from getdist import plots, MCSamples

conf = [0.683, 0.955, 0.997]
labels_xy = ['\\mathrm{X}', '\\mathrm{Y}']

# Build MCSamples for each star's XY
mcs_list = []
for xcol, ycol, sc in star_xy:
    mcs = MCSamples(samples=sc[[xcol, ycol]].values, names=[xcol, ycol], labels=labels_xy, ignore_rows=0.1)
    mcs.updateSettings({'contours': conf})
    mcs_list.append((xcol, ycol, mcs))

# Axis limits centred on star 1 position at chain index 1000
xcol0, ycol0, sc0 = star_xy[0]
ref_x = sc0[xcol0][1000]
ref_y = sc0[ycol0][1000]
lims = [ref_x - 25, ref_x + 25, ref_y - 25, ref_y + 25]

# Color sets (base, mid, inner) for secondary stars
secondary_color_sets = [
    ('darkgrey',  'silver',         'white'),
    ('darkblue',  'cornflowerblue', 'white'),
]

g = plots.getSubplotPlotter()
g.settings.num_plot_contours = 3
g.settings.shade_level_scale = 3
g.settings.fig_width_inch = 9
g.settings.linewidth = 0.5

# Star 1: default color, 3 contour levels
xcol, ycol, mcs = mcs_list[0]
g.plot_2d(mcs, xcol, ycol, filled=False, lims=lims)

# Secondary stars: layered contour triplet with distinct color sets
for idx, (xcol, ycol, mcs) in enumerate(mcs_list[1:]):
    base, mid, inner = secondary_color_sets[idx % len(secondary_color_sets)]
    g.settings.num_plot_contours = 3
    g.plot_2d(mcs, xcol, ycol, filled=False, colors=[base],  lims=lims)
    g.settings.num_plot_contours = 2
    g.plot_2d(mcs, xcol, ycol, filled=False, colors=[mid],   lims=lims)
    g.settings.num_plot_contours = 1
    g.plot_2d(mcs, xcol, ycol, filled=False, colors=[inner], lims=lims)

plt.text(ref_x - 22.28, ref_y + 20.58, 'MB11262', fontsize=35, fontweight='heavy', color='white')
plt.text(ref_x + 11.28, ref_y - 22.28, 'OSIRIS',  fontsize=35, fontweight='heavy', color='yellow')

plt.imshow(image_data, cmap='viridis', vmin=10, vmax=200, interpolation='nearest')

plt.xlabel('X [pix]', fontsize=28)
plt.ylabel('Y [pix]', fontsize=28)
plt.tick_params(which='both', length=12, width=2, labelsize=26, direction='in', right=True, top=True)

g.export(f'contour_plot_{star_suffix}.pdf')

# =============================================================================
# CORNER PLOT
# =============================================================================

burnin = int(0.1 * len(chains_all))

# Assemble corrected array: XY columns (already -1 corrected) + extra columns as-is
xy_parts   = [sc[[xcol, ycol]].values[burnin:] for xcol, ycol, sc in star_xy]
extra_part = [chains_extra.values[burnin:]] if extra_col_names else []
all_chains = np.column_stack(xy_parts + extra_part)

corner_labels = list(col_names)

# Read bestfit values (same column order and XY correction as all_chains)
bf_array = None
if os.path.exists(bf_fname):
    with open(bf_fname) as f:
        bf_col_names = f.readline().lstrip('#').split()
    bf_row = pd.read_table(bf_fname, sep=r'\s+', skiprows=1, names=bf_col_names, nrows=1)
    for xcol, ycol, _ in star_xy:
        bf_row[xcol] -= 1
        bf_row[ycol] -= 1
    xy_bf    = np.concatenate([bf_row[[xcol, ycol]].values[0] for xcol, ycol, _ in star_xy])
    extra_bf = bf_row[extra_col_names].values[0] if extra_col_names else np.array([])
    bf_array = np.concatenate([xy_bf, extra_bf])

# Drop constant columns
varying = np.std(all_chains, axis=0) > 0
all_chains    = all_chains[:, varying]
corner_labels = [l for l, v in zip(corner_labels, varying) if v]
bestfit_values = bf_array[varying] if bf_array is not None else None

with plt.rc_context({'font.size': 8, 'font.family': 'courier new'}):
    fig_corner = corner.corner(
        all_chains,
        labels=corner_labels,
        #quantiles=[0.1587, 0.5, 0.8413],
        truths=bestfit_values,
        show_titles=True,
        title_fmt='.3f',
        title_kwargs={'fontsize': 7, 'pad': 2},
        label_kwargs={'fontsize': 9},
        levels=conf,
        plot_contours=True,
        #plot_density=True,
        plot_datapoints=False,
        fill_contours=True,
        smooth=1.0,
        smooth1d=1.0,
        #color='steelblue',
        contourf_kwargs={'alpha': 0.001},
        hist_kwargs={'alpha': 0.99},
        max_n_ticks=3,
    )
    for ax in fig_corner.get_axes():
        ax.tick_params(labelsize=6, pad=2)
        ax.xaxis.label.set_size(9)
        ax.yaxis.label.set_size(9)
    #fig_corner.suptitle('MCMC Star Position Distributions', fontsize=10, y=1.01)
    fig_corner.savefig(f'corner_plot_{star_suffix}.pdf', bbox_inches='tight', dpi=150)
    plt.close(fig_corner)

# =============================================================================
# CHI^2 PIXEL PLOT
# =============================================================================

c_data = np.loadtxt(c_fname)

# Columns: X, Y, chi^2, intensity
X = c_data[:, 0]
Y = c_data[:, 1]
chi2 = c_data[:, 2]
intensity = c_data[:, 3]

# Create side-by-side plots
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# --- Left: chi^2 ---
sc1 = axes[0].scatter(X, Y, c=chi2, marker='s', s=50, edgecolor='None', cmap='viridis')
axes[0].set_title(r'$\chi^2$')
axes[0].set_xlabel('X')
axes[0].set_ylabel('Y')
plt.colorbar(sc1, ax=axes[0], label=r'$\chi^2$')

# --- Right: intensity ---
sc2 = axes[1].scatter(X, Y, c=intensity, marker='s', s=50, edgecolor='None', cmap='viridis')
axes[1].set_title('Intensity')
axes[1].set_xlabel('X')
axes[1].set_ylabel('Y')
plt.colorbar(sc2, ax=axes[1], label='Intensity')

for ax in axes:
    ax.set_xlim(lims[0], lims[1])
    ax.set_ylim(lims[2], lims[3])
    ax.ticklabel_format(style='plain')

plt.tight_layout()
plt.savefig(f'chi2_pixel_plots_{star_suffix}.pdf', dpi=150)
plt.close()
