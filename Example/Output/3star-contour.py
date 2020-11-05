import getdist
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image, ImageOps, ImageChops

img = Image.open("ds9.png")
img = ImageOps.flip(img)
img = ImageChops.offset(img,1,1)

plt.rcParams.update({'font.size': 110})
plt.rcParams["font.family"]="courier new"

colname = ['Separation', 'Flux_Ratio']
fname = 'new-mcmc_fit.dat'
chains1 = pd.read_table(fname, usecols=[0,1], sep='\s+', skiprows=1, names=colname)
chains2 = pd.read_table(fname, usecols=[2,3], sep='\s+', skiprows=1, names=colname)
chains3 = pd.read_table(fname, usecols=[4,5], sep='\s+', skiprows=1, names=colname)

from getdist import plots, MCSamples
labels = ['\\mathrm{X}', '\\mathrm{Y}']
samples1 = MCSamples(samples=chains1[['Separation', 'Flux_Ratio']].values, names=colname, labels=labels)
samples2 = MCSamples(samples=chains2[['Separation', 'Flux_Ratio']].values, names=colname, labels=labels)
samples3 = MCSamples(samples=chains3[['Separation', 'Flux_Ratio']].values, names=colname, labels=labels)

conf = [0.683, 0.955, 0.997]
samples1.updateSettings({'contours': conf})
samples2.updateSettings({'contours': conf})
samples3.updateSettings({'contours': conf})

g = plots.getSubplotPlotter()
g.settings.num_plot_contours = 3
g.settings.shade_level_scale = 3
g.settings.fig_width_inch = 8.5
g.settings.linewidth = 0.3
g.plot_2d(samples1,'Separation','Flux_Ratio', filled=False, lims=[610,640,360,390]) #30 x 30 pixels
g.plot_2d(samples2,'Separation','Flux_Ratio', filled=False, lims=[610,640,360,390]) #30 x 30 pixels
g.plot_2d(samples3,'Separation','Flux_Ratio', filled=False, colors=['skyblue'], lims=[610,640,360,390])
g.settings.num_plot_contours = 2
g.plot_2d(samples3,'Separation','Flux_Ratio', filled=False, colors=['lightblue'], lims=[610,640,360,390])
g.settings.num_plot_contours = 1
g.plot_2d(samples3,'Separation','Flux_Ratio', filled=False, colors=['white'], lims=[610,640,360,390])
plt.imshow(img)

g.export('c-overlay.pdf')
