# -*- coding: utf-8 -*-
"""
curve fitting tools
@author: Jussi Tiira
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
import seaborn as sns

sns.set_style('ticks')

GUNN_KINZER = (9.65, 10.30/9.65, 0.6)


def antidiagonal_identity(n):
    return np.matrix(np.identity(n))[::-1]


def antidiagonal_transpose(matrix):
    n = len(matrix)
    J = antidiagonal_identity(n)
    return J*matrix.T*J


class ClassProperty(property):
    def __get__(self, cls, owner):
        return self.fget.__get__(None, owner)()


class Fit:
    """parent for different fit types"""
    def __init__(self, x=None, y=None, x_unfiltered=None, y_unfiltered=None,
                 sigma=None, params=None, name='fit', xname='D',
                 flipped=False):
        self.params = params
        self._name = name
        self.x = x
        self.y = y
        self.x_unfiltered = x_unfiltered
        self.y_unfiltered = y_unfiltered
        self.fltr_upper_x = None
        self.fltr_upper_y = None
        self.fltr_lower_x = None
        self.fltr_lower_y = None
        self.sigma = sigma
        self.xname = xname
        self.flipped = flipped
        self.str_fmt = ''

    def __repr__(self):
        if self.params is None:
            paramstr = 'abcdefghij'[:self.str_fmt.count('%')]
        else:
            paramstr = ['{0:.3f}'.format(p) for p in self.params]
        s = self.str_fmt % tuple(paramstr)
        return s.replace('--', '+')

    def func(self, x, a=None):
        """Fit function. If no coefficients are given use stored ones."""
        if a is None:
            return self.func(x, *self.params)
        pass

    def penalty(self, params):
        """penalty function used by the cost function"""
        return 0

    def plot(self, xmin=None, xmax=None, samples=1000, ax=None, label=None,
             linewidth=2, source_style=None, unfiltered=False, hide_filter=False,
             source_kwargs={}, **kwargs):
        """Plot fit curve and fitted data."""
        if ax is None:
            ax = plt.gca()
        if self.params is None:
            return ax
        if xmax is None:
            if self.x is None or len(self.x)<1:
                xmax = 10
            else:
                xmax = self.x.max()
        if xmin is None:
            if self.x is None or len(self.x)<1:
                xmin = 0
            else:
                xmin = self.x.min()
        x = np.linspace(xmin, xmax, samples)
        y = [self.func(xi, *self.params) for xi in x]
        if label is None:
            label = r'$' + str(self) + r'$'
        ax.plot(x, y, label=label, linewidth=linewidth, **kwargs)
        if unfiltered:
            x = self.x_unfiltered
            y = self.y_unfiltered
            if self.fltr_upper_x is not None and not hide_filter:
                fltr_kws = {'where': 'post',
                            'linestyle': 'dotted',
                            'color': 'red'}
                ax.step(self.fltr_upper_x, self.fltr_upper_y, **fltr_kws)
                ax.step(self.fltr_lower_x, self.fltr_lower_y, **fltr_kws)
        else:
            x = self.x
            y = self.y
        if source_style=='raw':
            ax.scatter(x, y, **source_kwargs)
        elif source_style=='kde':
            sns.kdeplot(x, y, ax=ax, shade=True, shade_lowest=False,
                        bw=.01, **source_kwargs)
        elif source_style=='hex':
            ax.hexbin(x, y, **source_kwargs)
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
        if self.x is None or self.y is None:
            return
        if self.sigma is not None:
            kwargs['sigma'] = self.sigma
        elif self.flipped:
            x = self.y
            y = self.x
        else:
            x = self.x
            y = self.y
        params, cov = curve_fit(self.func, x, y, **kwargs)
        if self.flipped:
            params = self.flip_params(params)
            cov = antidiagonal_transpose(cov)
        if store_params:
            self.params = params
            self.cov = cov
        return params, cov

    def flip_params(self, params):
        return params

    def perr(self, cov=None):
        """standard errors of x and y"""
        if cov is None:
            cov = self.cov
        return np.sqrt(np.diag(cov))

    def is_good(self):
        return True

    @ClassProperty
    @classmethod
    def name(cls):
        return cls.__name__.lower()


class LinFit(Fit):
    def __init__(self, params=None, **kwargs):
        super().__init__(params=params, name='linfit', **kwargs)
        self.str_fmt = '%s' + self.xname + '+%s'

    def func(self, x, a=None, b=None):
        if a is None:
            return self.func(x, *self.params)
        return a*x+b

    def find_fit(self, store_params=True, **kwargs):
        if self.x is None or self.y is None:
            return
        x = self.x
        y = self.y
        a, b, r_value, p_value, std_err = linregress(x, y)
        self.params = (a, b)
        return self.params


class ExpFit(Fit):
    """exponential fit of form a*(1-b*exp(-c*D))"""
    def __init__(self, params=None, **kwargs):
        super().__init__(params=params, name='expfit', **kwargs)
        self.quess = (1., 1., 1.)
        self.str_fmt = '%s(1-%s\exp(-%s' + self.xname + '))'

    def func(self, x, a=None, b=None, c=None):
        if a is None:
            return self.func(x, *self.params)
        if self.flipped:
            return -1/c*np.log((1-x/a)/b)
        return a*(1-b*np.exp(-c*x))

    def penalty(self, params):
        return 0
        return 1000*(max(0, 0.1-params[1]) + max(0, 0.4-params[2]))


class LogFit(Fit):
    def __init__(self, params=None, **kwargs):
        super().__init__(params=params, name='logfit', **kwargs)

    def func(self, x, a=None, b=None):
        if a is None:
            return self.func(x, *self.params)
        return a*np.log10(x) + b


class ExponentialFit(Fit):
    def __init__(self, params=None, base=10, **kwargs):
        super().__init__(params=params, name='logfit', **kwargs)
        self.base = base
        self.str_fmt = '%s' + str(base) + '^{%s' + self.xname + '}'

    def func(self, x, a=None, b=None):
        if a is None:
            return self.func(x, *self.params)
        return a*self.base**(b*x)


class PolFit(Fit):
    """power law fit of form a*D**b"""
    def __init__(self, params=None, **kwargs):
        super().__init__(params=params, name='polfit', **kwargs)
        self.quess = (1., 1.)
        self.str_fmt = '%s' + self.xname + '^{%s}'

    def func(self, x, a=None, b=None):
        if a is None:
            return self.func(x, *self.params)
        return a*x**b

    def penalty(self, params):
        #return 0
        return 1000*max(0, 0.2-params[1])

    def flip_params(self, params):
        return np.array([params[0]**(-1/params[1]), 1/params[1]])

    def is_good(self):
        a = self.params[0]
        b = self.params[1]
        if (a < 0) or (b < 0) or (b > 1):
            return False
        return True

gunn_kinzer = ExpFit(params=GUNN_KINZER)
