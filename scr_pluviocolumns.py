# -*- coding: utf-8 -*-
"""
@author: jussitii
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt

e = EventsCollection('cases/cases_of_interest.csv', '%d.%m. %H:%M')
e.autoimport_data(autoshift=False, rule='6min')

p4 = e.events.pluvio400[0].pluvio
p4.shift_periods = -5
p4.acc('1min').tshift(periods=5).plot(label='Bucket NRT filtered')
raw = p4.data.bucket_nrt - p4.data.bucket_nrt[0]
rawrt = p4.data.bucket_rt - p4.data.bucket_rt[0]
rawacc = p4.data.acc_tot_nrt - p4.data.acc_tot_nrt[0]
rawacc = p4.data.acc_rt.cumsum()
raw.tshift(periods=-5,freq='1min').plot(label='Bucket NRT -5min')
rawrt.plot(label='Bucket RT')
rawacc.plot(label='Accumulation RT/NRT')
rawacc.plot(label='Accumulated total NRT')
plt.legend(loc='lower right')