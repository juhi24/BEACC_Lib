import pandas as pd
import matplotlib.pyplot as plt
import read
from glob import glob

pluvio_files = glob('../DATA/Pluvio200/*.txt')
dsd_files = glob('../DATA/PIP/a_DSD_Tables/00420140221*.dat')
pipv_files = glob('../DATA/PIP/a_Velocity_Tables/00420140221/*2.dat')

pluvio = read.Pluvio(pluvio_files)
dsd = read.PipDSD(dsd_files)
pipv = read.PipV(pipv_files)

f1 = plt.figure(figsize=(6,5), dpi=100)

pipv.data.between_time('23:30','00:00').plot(x='Wad_Dia',y='vel_v',style='+')
plt.savefig('../vel.png',bbox_inches='tight')
plt.clf()

pipv.data.between_time('23:15','23:20').vel_v.hist(bins=100)
plt.savefig('../vel_hist.png')
plt.clf()
