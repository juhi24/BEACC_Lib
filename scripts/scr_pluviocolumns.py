# -*- coding: utf-8 -*-
"""
@author: jussitii
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt

e = EventsCollection('cases/cases_of_interest.csv', '%d.%m. %H:%M')
e.autoimport_data(autoshift=False, rule='6min')

plt.close('all')

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    plt.figure()
    p4 = c.pluvio
    tshift = -0
    #p4.shift_periods = 0
    acc_shift = 0.00
    #p4.acc('1min').tshift(periods=5).plot(label='Bucket NRT filtered')
    raw = p4.data.bucket_nrt - p4.data.bucket_rt[0]
    rawrt = p4.data.bucket_rt - p4.data.bucket_rt[0] + acc_shift
    rawacc = p4.data.acc_tot_nrt - p4.data.acc_tot_nrt[0]
    rawacc = p4.data.acc_rt.cumsum()
    sign = ''
    if tshift>-0.1:
        sign = '+'
    raw.tshift(periods=tshift,freq='1min').plot(label='Bucket NRT ' + sign + str(tshift) + 'min')
    rawrt.plot(label='Bucket RT')
    #rawacc.plot(label='Accumulation RT/NRT')
    #rawacc.plot(label='Accumulated total NRT')
    plt.title(p4.name + ' ' + str(raw.index[0].date()))
    plt.legend(loc='upper left')