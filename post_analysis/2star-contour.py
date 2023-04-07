import getdist
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image, ImageOps, ImageChops

import pdb

img = Image.open("ds9.png")
img = ImageOps.flip(img)
img = ImageChops.offset(img,1,1)

plt.rcParams.update({'font.size': 110})
plt.rcParams["font.family"]="courier new"

colname1 = ['X1_CENTER', 'Y1_CENTER']
colname2 = ['X2_CENTER', 'Y2_CENTER']
fname = 'accepted.mcmc'
chains1 = pd.read_table(fname, usecols=[0,1], sep='\s+', skiprows=1, names=colname1)
chains2 = pd.read_table(fname, usecols=[2,3], sep='\s+', skiprows=1, names=colname2)

from getdist import plots, MCSamples
labels = ['\\mathrm{X}', '\\mathrm{Y}']
samples1 = MCSamples(samples=chains1[['X1_CENTER', 'Y1_CENTER']].values, names=colname1, labels=labels, ignore_rows=0.1) 
samples2 = MCSamples(samples=chains2[['X2_CENTER', 'Y2_CENTER']].values, names=colname2, labels=labels, ignore_rows=0.1)

conf = [0.683, 0.955, 0.997]
samples1.updateSettings({'contours': conf})
samples2.updateSettings({'contours': conf})

g = plots.getSubplotPlotter()
g.settings.num_plot_contours = 3
g.settings.shade_level_scale = 3
g.settings.fig_width_inch = 8.5
g.settings.linewidth = 0.1
g.plot_2d(samples1,'X1_CENTER','Y1_CENTER', filled=False, lims=[chains1['X1_CENTER'][1000]-20, chains1['X1_CENTER'][1000]+20,chains1['Y1_CENTER'][1000]-20,chains1['Y1_CENTER'][1000]+20])
g.plot_2d(samples2,'X2_CENTER','Y2_CENTER', filled=False, colors=['skyblue'], lims=[chains1['X1_CENTER'][1000]-20, chains1['X1_CENTER'][1000]+20,chains1['Y1_CENTER'][1000]-20,chains1['Y1_CENTER'][1000]+20])
g.settings.num_plot_contours = 2
g.plot_2d(samples2,'X2_CENTER','Y2_CENTER', filled=False, colors=['lightblue'], lims=[chains1['X1_CENTER'][1000]-20, chains1['X1_CENTER'][1000]+20,chains1['Y1_CENTER'][1000]-20,chains1['Y1_CENTER'][1000]+20])
g.settings.num_plot_contours = 1
g.plot_2d(samples2,'X2_CENTER','Y2_CENTER', filled=False, colors=['white'], lims=[chains1['X1_CENTER'][1000]-20, chains1['X1_CENTER'][1000]+20,chains1['Y1_CENTER'][1000]-20,chains1['Y1_CENTER'][1000]+20])
plt.imshow(img)

g.export('c-overlay.pdf')
