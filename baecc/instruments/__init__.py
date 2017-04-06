# coding: utf-8
from __future__ import absolute_import, division, print_function
import gc
import baecc
from baecc.instruments.base import InstrumentData, PrecipMeasurer, datafilelistloop
from baecc.instruments import pluvio
from baecc.instruments import pip_psd
from baecc.instruments import pip_v
from baecc.instruments import radar

# PIP-observed particle size to the volume equivalent diameter
PHI = 0.92


def batch_import(dtstrlist, datadir=baecc.DATA_DIR, use_radar=False):
    """Read ASCII data according to a datestring pattern."""
    pipv_files = datafilelistloop(pip_v.SUBPATH, dtstrlist, datadir=datadir)
    dsd_files = datafilelistloop(pip_psd.SUBPATH, dtstrlist, datadir=datadir)
    pluvio200_files = datafilelistloop(pluvio.P200_SUBPATH, dtstrlist, datadir=datadir)
    pluvio400_files = datafilelistloop(pluvio.P400_SUBPATH, dtstrlist, datadir=datadir)
    if use_radar:
        xsacr_files = datafilelistloop(radar.SUBPATH, [('XSACR', 'xsacr', dtstr) for dtstr in dtstrlist],
                                   datadir=datadir)
        kasacr_files = datafilelistloop(radar.SUBPATH, [('KASACR', 'kasacr', dtstr) for dtstr in dtstrlist],
                                    datadir=datadir)
        kazr_files = datafilelistloop(radar.SUBPATH, [('KAZR', 'kazrge', dtstr) for dtstr in dtstrlist],
                                  datadir=datadir)
        mwacr_files = datafilelistloop(radar.SUBPATH, [('MWACR', 'mwacr', dtstr) for dtstr in dtstrlist],
                                   datadir=datadir)
    pluvio200 = pluvio.Pluvio(pluvio200_files)
    pluvio400 = pluvio.Pluvio(pluvio400_files)
    pipv = pip_v.PipV(pipv_files)
    dsd = pip_psd.PipPSD(dsd_files)
    if use_radar:
        xsacr = radar.Radar(xsacr_files)
        kasacr = radar.Radar(kasacr_files)
        kazr = radar.Radar(kazr_files)
        mwacr = radar.Radar(mwacr_files)
        return {'vel': pipv, 'dsd': dsd, 'pluvio200': pluvio200,
                'pluvio400': pluvio400, 'xsacr': xsacr, 'kasacr': kasacr,
                'kazr': kazr, 'mwacr': mwacr}
    return {'vel': pipv, 'dsd': dsd, 'pluvio200': pluvio200,
            'pluvio400': pluvio400}


def batch_create_hdf(instrdict=None, datadir=baecc.DATA_DIR,
                     hdf_file=baecc.H5_PATH, dtstrlist=('20140[2-3]??')):
    """Read ASCII data and export to hdf5."""
    if instrdict is None:
        instrdict = {pip_v.SUBPATH: pip_v.PipV,
                     pip_psd.SUBPATH: pip_psd.PipPSD,
                     pluvio.P200_SUBPATH: pluvio.Pluvio,
                     pluvio.P400_SUBPATH: pluvio.Pluvio}
    for key in instrdict:
        instr = instrdict[key].from_raw(dtstrlist, subpath=key, datadir=datadir)
        instr.to_hdf(filename=hdf_file)
        del(instr)
        gc.collect()
