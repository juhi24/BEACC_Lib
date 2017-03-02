# coding: utf-8
import numpy as np
import pandas as pd
from baecc.instruments import InstrumentData
from pytmatrix import psd

RHO_0 = 1.225 # sea level air density, from http://www.aviation.ch/tools-atmosphere.asp
RHO_W = 1000.0


def ar(d):
    return 1/aux.dsr_thurai_2007(d)

def drop_zeros_rows(df):
    return df[(df.T != 0).any()]

class PSD(InstrumentData):
    """General PSD"""
    # TODO: fix overlap with PipDSD
    # TODO: make PipPSD inherit this
    def __init__(self, data=None, binwidth=None, bin_v=None):
        if data is not None:
            self.data = data.resample('1min').fillna(0)
        self.binwidth = binwidth
        self.bin_v = pd.Series(bin_v, index=data.columns.values)
        self.drop_empty = True

    def plot(self, **kwargs):
        """wrapper for plot_psd"""
        return plot_psd(self.good_data(), **kwargs)

    def binwidth_df(self):
        return pd.DataFrame(data=self.binwidth, index=self.bin_cen()).T

    def good_data(self, drop_empty=True):
        data = self.data
        if self.drop_empty and drop_empty:
            data = drop_zeros_rows(data)
        return data

    def intensity(self):
        return self.sum_over_d(self.rr)

    def v(self, d):
        return self.bin_v[self.bin_select(d)]

    def n(self, d):
        return self.good_data()[d]

    def bin_cen(self):
        return self.good_data().columns.values

    def bin_lower(self):
        return self.bin_cen() - 0.5*self.binwidth

    def bin_upper(self):
        return self.bin_lower() + self.binwidth

    def bin_edge(self):
        return np.append(self.bin_lower(), self.bin_upper()[-1])

    def bin_edge_df(self):
        edge = pd.DataFrame(data=[self.bin_lower(), self.bin_upper()]).T
        edge.columns = ['lower', 'upper']
        edge.index = self.bin_cen()
        return edge

    def bin_select(self, d):
        for i, edge in self.bin_edge_df().iterrows():
            if d > edge.lower and d < edge.upper:
                return edge.name
        return

    def series_zeros(self):
        s = self.good_data()[self.good_data().columns[0]]*0
        s.name = 'series'
        return s

    def sum_over_d(self, func, **kwargs):
        dD = self.binwidth
        result = self.series_zeros()
        for d in self.bin_cen():
            result = result.add(func(d, **kwargs)*dD[d], fill_value=0)
        return result

    def rr(self, d):
        return 3.6e-3*np.pi/6*(ar(d)*d)**2*1/ar(d)*d*self.v(d)*self.n(d)

    def to_tm(self, data=None):
        if data is None:
            data = self.good_data().mean()
        return psd.BinnedPSD(bin_edges=self.bin_edge(),
                             bin_psd=data.values)

    def to_tm_series(self, resample=None):
        if resample is None:
            data = self.good_data()
        else:
            data = self.good_data(drop_empty=False).resample(resample, how=np.mean,
                                                             closed='right',
                                                             label='right')
        return data.apply(self.to_tm, axis=1)
