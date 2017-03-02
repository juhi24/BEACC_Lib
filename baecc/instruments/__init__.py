# coding: utf-8
import gc
from baebb import H5_PATH, DATA_DIR
from baecc.instruments.pluvio import Pluvio
from baecc.instruments.pip_psd import PipPSD
from baecc.instruments.pip_v import PipV
from baecc.instruments.radar import Radar
from baecc.instruments.base import InstrumentData, PrecipMeasurer
from baecc.instruments.tools import datafilelistloop

PIPV_SUBPATH = 'PIP/a_Velocity_Tables/004%s/*2.dat'
PSD_SUBPATH = 'PIP/a_DSD_Tables/004%s*.dat'
P200_SUBPATH = 'Pluvio200/pluvio200_??_%s*.txt'
P400_SUBPATH = 'Pluvio400/pluvio400_??_%s*.txt'
RADAR_SUBPATH = 'Radar/%s/tmp%s*M1.a1.%s.*'

# PIP-observed particle size to the volume equivalent diameter
PHI = 0.92


def batch_import(dtstrlist, datadir=DATA_DIR, radar=False):
    """Read ASCII data according to a datestring pattern."""
    pipv_files = datafilelistloop(PIPV_SUBPATH, dtstrlist, datadir=datadir)
    dsd_files = datafilelistloop(PSD_SUBPATH, dtstrlist, datadir=datadir)
    pluvio200_files = datafilelistloop(P200_SUBPATH, dtstrlist, datadir=datadir)
    pluvio400_files = datafilelistloop(P400_SUBPATH, dtstrlist, datadir=datadir)
    if radar:
        xsacr_files = datafilelistloop(RADAR_SUBPATH, [('XSACR', 'xsacr', dtstr) for dtstr in dtstrlist],
                                   datadir=datadir)
        kasacr_files = datafilelistloop(RADAR_SUBPATH, [('KASACR', 'kasacr', dtstr) for dtstr in dtstrlist],
                                    datadir=datadir)
        kazr_files = datafilelistloop(RADAR_SUBPATH, [('KAZR', 'kazrge', dtstr) for dtstr in dtstrlist],
                                  datadir=datadir)
        mwacr_files = datafilelistloop(RADAR_SUBPATH, [('MWACR', 'mwacr', dtstr) for dtstr in dtstrlist],
                                   datadir=datadir)
    pluvio200 = Pluvio(pluvio200_files)
    pluvio400 = Pluvio(pluvio400_files)
    pipv = PipV(pipv_files)
    dsd = PipPSD(dsd_files)
    if radar:
        xsacr = Radar(xsacr_files)
        kasacr = Radar(kasacr_files)
        kazr = Radar(kazr_files)
        mwacr = Radar(mwacr_files)
        return {'vel': pipv, 'dsd': dsd, 'pluvio200': pluvio200,
                'pluvio400': pluvio400, 'xsacr': xsacr, 'kasacr': kasacr,
                'kazr': kazr, 'mwacr': mwacr}
    return {'vel': pipv, 'dsd': dsd, 'pluvio200': pluvio200,
            'pluvio400': pluvio400}


def batch_create_hdf(instrdict=None, datadir=DATA_DIR, hdf_file=H5_PATH,
                     dtstrlist=('20140[2-3]??')):
    """Read ASCII data and export to hdf5."""
    if instrdict is None:
        instrdict = {PIPV_SUBPATH: PipV,
                     PSD_SUBPATH: PipPSD,
                     P200_SUBPATH: Pluvio,
                     P400_SUBPATH: Pluvio}
    for key in instrdict:
        instr = instrdict[key].from_raw(dtstrlist, subpath=key, datadir=datadir)
        instr.to_hdf(filename=hdf_file)
        del(instr)
        gc.collect()
