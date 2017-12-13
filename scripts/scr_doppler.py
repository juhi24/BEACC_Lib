# coding: utf-8
import pandas as pd
import matplotlib.pyplot as plt
import read
from glob import glob
from datetime import datetime, timedelta

#pluvio_files = glob('../DATA/Pluvio200/*.txt')
#dsd_files = glob('../DATA/PIP/a_DSD_Tables/00420140221*.dat')
pipv_files = glob('../DATA/PIP/a_Velocity_Tables/00420140221/*2.dat')

#pluvio = read.Pluvio(pluvio_files)
#dsd = read.PipDSD(dsd_files)
pipv = read.PipV(pipv_files)

f1 = plt.figure(figsize=(6,5), dpi=100)

dt = datetime(2014,2,21,22,45)
for i in range(12):
   dt_next = dt + timedelta(minutes=5)
   timestr = str(dt.time())
   nextstr = str(dt_next.time())
   hh = timestr[0:2]
   mm = timestr[3:5]

   pipv.data.between_time(timestr,nextstr).plot(x='Wad_Dia',y='vel_v',style='+')
   partcount = pipv.data.between_time(timestr,nextstr).Part_ID.count()
   plt.axis([0,12,0,5])
   plt.title('PIP: particle size vs. fall velocity %s - %s' % (timestr[0:5],nextstr[0:5]))
   plt.text(0.5, 4.6, 'particle count: %s' % str(partcount))
   plt.xlabel('D (mm)')
   plt.ylabel('Vertical velocity (m/s)')
   plt.savefig('../vel%s%s.png' % (hh,mm))
   plt.clf()

   pipv.data.between_time(timestr,nextstr).vel_v.hist(bins=100)
   plt.axis([0.5,4,0,450])
   plt.title('PIP: fall velocity histogram %s - %s' % (timestr[0:5],nextstr[0:5]))
   plt.xlabel('Vertical velocity (m/s)')
   plt.ylabel('Particle count')
   plt.savefig('../vel_hist%s%s.png' % (hh,mm))
   plt.clf()

   dt = dt_next
