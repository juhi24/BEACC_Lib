# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
import snowfall as sf
import numpy as np
from scipy.special import gamma

def rho(alpha, beta, d0, b=0.2):
    return 6/np.pi*alpha*gamma(beta+b+1)/gamma(b+4)*(d0/3.67)**(beta-3)