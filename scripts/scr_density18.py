# coding: utf-8
from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import pandas as pd
from os import path
from glob import glob
from baecc.instruments import pluvio, pip_psd, pip_v
from j24 import home, ensure_join


RESULTS_DIR = ensure_join(home(), 'results', 'density18')


def make_hdf():
    datadir = path.join(home(), 'DATA', 'HYY_DenData')
    p200dir = path.join(datadir, 'Pluvio200')
    p400dir = path.join(datadir, 'Pluvio400')
    pipdir = path.join(datadir, 'PIP')
    psddir = path.join(pipdir, 'f_1_4_DSD_Tables_ascii')
    vdir = path.join(pipdir, 'f_2_2_Velocity_Tables')
    h5file = path.join(datadir, 'annakaisa15-18.h5')
    pluv_paths = dict(pluvio200=p200dir, pluvio400=p400dir)
    instr = dict()
    # read pluvio data
    for p in ['pluvio200', 'pluvio400']:
        fnames = glob(path.join(pluv_paths[p], '*.txt'))
        fnames.sort()
        pluv = pluvio.Pluvio(fnames, name=p)
        selection = pluv.data.i_rt.apply(lambda x: type(x)==float)
        pluv.data = pluv.data[selection].astype(float)
        instr[p] = pluv
    # read pip data
    vfiles = glob(path.join(vdir, '*', '*.dat'))
    vfiles.sort()
    psdfiles = glob(path.join(psddir, '*.dat'))
    psdfiles.sort()
    instr['vel'] = pip_v.PipV(vfiles)
    instr['dsd'] = pip_psd.PipPSD(psdfiles)
    for i in instr.values():
        i.to_hdf(filename=h5file)


def save_density(rho, savedir=RESULTS_DIR):
    filename = rho.index[0].date().isoformat() + '.csv'
    filepath = path.join(savedir, filename)
    rho.to_csv(filepath)


if __name__ == '__main__':
    casesdir = path.join(home(), '.baecc', 'cases')
    cfile = path.join(casesdir, 'annakaisa15-18.csv')
    dates = pd.read_csv(cfile, parse_dates=[0,1])
    cases = []
    rhos = []
    for i, row in dates.iterrows():
        c = c2.between_datetime(row.start, row.end)
        c.use_cache = False
        try:
            rhos.append(c.density())
        except Exception:
            pass
        cases.append(c)
    for rho in rhos:
        save_density(rho)

