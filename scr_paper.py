# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
import gc
runfile('scr_case_overview.py')
gc.collect()
runfile('scr_histograms.py')
gc.collect()
runfile('scr_psd_plots.py')
gc.collect()
runfile('scr_d0-rho_combined.py')
gc.collect()
runfile('scr_vfits_density_ranges_combined.py')
gc.collect()