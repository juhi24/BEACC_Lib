# coding: utf-8
import numpy as np
import pandas as pd
import datetime
from baecc import instruments


class PipParticles(instruments.InstrumentData):
    """PIP particle tables"""
    def __init__(self, filenames=None, dt_start=None, dt_end=None, **kwargs):
        InstrumentData.__init__(self, filenames, **kwargs)
        self.name = 'pip_part'
        dtype = {'Year': np.int32, 'Month': np.int32, 'Day': np.int32,
                 'Hr': np.int32, 'Min': np.int32, 'Sec': np.int32}
        if self.data.empty:
            print('Reading PIP particle data...')
            for filename in filenames:
                newdata = pd.read_csv(filename, delim_whitespace=True,
                                      skiprows=[0, 1, 2, 3, 4, 5, 6, 7, 9],
                                      index_col='datetime',
                                      parse_dates={'datetime': ['Year',
                                                                'Month',
                                                                'Day', 'Hr',
                                                                'Min', 'Sec']},
                                      date_parser=datetime.datetime,
                                      dtype=dtype)
                self.data = self.data.append(newdata)
            self.finish_init(dt_start, dt_end)
