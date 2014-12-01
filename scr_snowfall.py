# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt

e = EventsCollection('cases/cases_of_interest.csv', '%d.%m. %H:%M')
e.autoimport_data(autoshift=False, autobias=False, rule='6min')

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.pluvio.shift_periods = -6
    #c.minimize_lsq()
    #c.plot(pip=False)
    #c.pipv.plots(save=True)
    
c.ab=(0.1,2.1)
vels = c.pipv.v(0.375, rule=c.pluvio.grouper())

p=e.events.pluvio200[2].pluvio
d=p.data
amount=d.acc_nrt[d.acc_nrt>0]
a=amount.between_time('22:13', '23:30')

# delta method
v=e.events.pluvio200[2].pipv
vel = v.good_data().vel_v
tdelta = pd.Series(a.index.to_pydatetime(), index=a.index).diff()
tdelta[0] = pd.datetools.timedelta(minutes=3)

# groupby method
am = d.acc_nrt.between_time('22:13', '23:30')
group = am.astype(bool).astype(int).cumsum().shift(1).fillna(0)
group.name = 'group'
gdf=pd.DataFrame(group)
vdf=pd.DataFrame(vel)
vgrouped=pd.merge(vdf,gdf,left_index=True,right_index=True)