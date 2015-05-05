# -*- coding: utf-8 -*-
"""
curve fitting tools
@author: Jussi Tiira
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

GUNN_KINZER = (9.65, 10.30/9.65, 0.6)

class Fit:
    """parent for different fit types"""
    def __init__(self, x=None, y=None, sigma=None, params=None, name='fit',
                 xname='D'):
        self.params = params
        self.name = name
        self.x = x
        self.y = y
        self.sigma = sigma
        self.xname = xname

    def func(self, x, a=None):
        """Fit function. If no coefficients are given use stored ones."""
        if a is None:
            return self.func(x, *self.params)
        pass

    def penalty(self, params):
        """penalty function used by the cost function"""
        return 0

    def plot(self, xmax=None, samples=1000, ax=None, label=None,
             source_data=False, marker='ro', linewidth=2, **kwargs):
        """Plot fit curve and fitted data."""
        if ax is None:
            ax = plt.gca()
        if self.params is None:
            return ax
        if xmax is None:
            if self.x is None:
                xmax = 10
            else:
                xmax = x.max()
        x = np.linspace(0, xmax, samples)
        y = [self.func(xi, *self.params) for xi in x]
        if label is None:
            label = r'$' + str(self) + r'$'
        ax.plot(x, y, label=label, linewidth=linewidth, **kwargs)
        if source_data:
            ax.plot(self.x, self.y, marker)
        return ax

    def cost(self, params, xarr, yarr, sigarr):
        """Cost function that can be used to find fit coefs by minimization."""
        cost = 0
        for i, x in enumerate(xarr):
            y = yarr[i]
            sig = sigarr[i]
            cost += 1/sig**2*(y - self.func(x, *params))**2 + self.penalty(params)
        return cost

    def find_fit(self, store_params=True, **kwargs):
        selection = pd.DataFrame([self.x.notnull(), self.y.notnull()]).all()
        x = self.x[selection]
        y = self.y[selection]
        if self.sigma is not None:
            kwargs['sigma'] = self.sigma[selection]
        params, cov = curve_fit(self.func, x, y, **kwargs)
        if store_params:
            self.params = params
        return params, cov

class ExpFit(Fit):
    """exponential fit of form a*(1-b*exp(-c*D))"""
    def __init__(self, params=None, **kwargs):
        super().__init__(params=params, name='expfit', **kwargs)
        self.quess = (1., 1., 1.)

    def __repr__(self):
        if self.params is None:
            paramstr = 'abc'
        else:
            paramstr = ['{0:.3f}'.format(p) for p in self.params]
        s = '%s(1-%s\exp(-%s%s))' % (paramstr[0], paramstr[1], paramstr[2],
                                     self.xname)
        return s.replace('--', '+')

    def func(self, x, a=None, b=None, c=None):
        if a is None:
            return self.func(x, *self.params)
        return a*(1-b*np.exp(-c*x))

    def penalty(self, params):
        return 0
        return 1000*(max(0, 0.1-params[1]) + max(0, 0.4-params[2]))

class PolFit(Fit):
    """polynomial fit of form a*D**b"""
    def __init__(self, params=None, **kwargs):
        super().__init__(params=params, name='polfit', **kwargs)
        self.quess = (1., 1.)

    def __repr__(self):
        if self.params is None:
            paramstr = 'ab'
        else:
            paramstr = ['{0:.3f}'.format(p) for p in self.params]
        return '%s%s^{%s}' % (paramstr[0], self.xname, paramstr[1])

    def func(self, x, a=None, b=None):
        if a is None:
            return self.func(x, *self.params)
        return a*x**b

    def penalty(self, params):
        #return 0
        return 1000*max(0, 0.2-params[1])

gunn_kinzer = ExpFit(params=GUNN_KINZER)
