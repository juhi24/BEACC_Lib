# -*- coding: utf-8 -*-
"""
curve fitting tools
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.stats import linregress
import seaborn as sns

GUNN_KINZER = (9.65, 10.30/9.65, 0.6)


def exp10(x):
    return int(np.floor(np.log10(abs(x))))


def num2tex(x, sig=3, sci_thres=2):
    exp = exp10(x)
    if abs(exp) > sci_thres:
        return num2tex_e(x, sig=sig)
    return '{0:.{sig}f}'.format(x, sig=sig)


def num2tex_e(x, sig=3):
    """return number as latex string in scientific format"""
    return '{0:.{sig}f} \mathrm{{e}}{{{1:-}}}'.format(*frexp10(x), sig=sig)


def num2tex_sci(x, sig=3):
    """return number as latex string in scientific format"""
    return '{0:.{sig}f} \\times 10^{{{1}}}'.format(*frexp10(x), sig=sig)


def frexp10(x):
    exp = exp10(x)
    return x / 10**exp, exp


def set_plot_style(tickdirection='in', **kws):
    styledict = {'xtick.direction': tickdirection,
                 'ytick.direction': tickdirection}
    styledict.update(kws)
    sns.set_style('ticks', styledict)


def antidiagonal_identity(n):
    return np.matrix(np.identity(n))[::-1]


def antidiagonal_transpose(matrix):
    n = len(matrix)
    J = antidiagonal_identity(n)
    return J*matrix.T*J


def logn(x, n=10):
    """n base logarithm using base change rule"""
    return np.log(n**x)/np.log(n)


class ClassProperty(property):
    def __get__(self, cls, owner):
        return self.fget.__get__(None, owner)()


class Fit:
    """parent for different fit types"""
    def __init__(self, x=None, y=None, x_unfiltered=None, y_unfiltered=None,
                 sigma=None, params=None, name='fit', xname='D',
                 flipped=False, disp_scale=[], use_latex_fmt=False):
        self.params = params
        self._name = name
        self.x = x
        self.y = y
        if x_unfiltered is None:
            self.x_unfiltered = x
            self.y_unfiltered = y
        else:
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
        self.disp_scale = np.array(disp_scale)
        self.use_latex_fmt = use_latex_fmt

    def __repr__(self):
        if self.params is None:
            paramstr = 'abcdefghijklmnopqrstuvwxyz'[:self.str_fmt.count('%')]
        else:
            params = np.array(self.params)
            if self.disp_scale.size == params.size:
                params = params*self.disp_scale
            try:
                if self.use_latex_fmt:
                    paramstr = [num2tex(p) for p in params]
                else:
                    paramstr = ['{0:.3f}'.format(p) for p in params]
            except AttributeError:
                paramstr = ['{0:.3f}'.format(p) for p in params]
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
             source_kws={}, **kwargs):
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
            ax.scatter(x, y, **source_kws)
        elif source_style=='kde':
            sns.kdeplot(x, y, ax=ax, shade=True, shade_lowest=False,
                        bw=.01, **source_kws)
        elif source_style=='hex':
            ax.hexbin(x, y, **source_kws)
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
        if self.flipped:
            x = self.y
            y = self.x
        else:
            x = self.x
            y = self.y
        params, cov = optimize.curve_fit(self.func, x, y, **kwargs)
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

    def xy(self, filtered=True):
        if filtered:
            return self.x, self.y
        else:
            return self.x_unfiltered, self.y_unfiltered

    def residuals(self, filtered=True, **kws):
        x, y = self.xy(filtered=filtered)
        return y - self.func(x, **kws)

    def sstot(self, filtered=True):
        """total sum of squares"""
        x, y = self.xy(filtered=filtered)
        return sum((y-y.mean())**2)

    def ssres(self, filtered=True, **kws):
        """residual sum of squares"""
        x, y = self.xy(filtered=filtered)
        return sum((y-self.func(x, **kws))**2)

    def rsq(self, filtered=True, **kws):
        """coefficient of determination, R**2"""
        return 1-self.ssres(filtered=filtered, **kws)/self.sstot(filtered=filtered)

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
        a, b, self.rvalue, self.pvalue, self.stderr = linregress(x, y)
        if store_params:
            self.params = (a, b)
        return self.params


class ExpFit(Fit):
    """exponential fit of form a*(1-b*exp(-c*x))"""
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
        self.str_fmt = '%s \\times' + str(base) + '^{%s ' + self.xname + '}'

    def func(self, x, a=None, b=None):
        if a is None:
            return self.func(x, *self.params)
        return a*self.base**(b*x)

    def find_fit(self, store_params=True, linreg=True, **kwargs):
        if not linreg:
            return super().find_fit(store_params=store_params, **kwargs)
        b, loga, rvalue, pvalue, stderr = self.find_fit_linregress()
        a = self.base**loga
        if store_params:
            self.params = (a, b)
            self.rvalue = rvalue
            self.pvalue = pvalue
            self.stderr = stderr
        return (a, b)

    def find_fit_linregress(self):
        basen = False
        if self.base == np.e:
            logfunc = np.log
        elif self.base == 10:
            logfunc = np.log10
        elif self.base == 2:
            logfunc = np.log2
        else:
            basen = True
        if basen:
            logy = logn(self.y, n=self.base)
        else:
            logy = logfunc(self.y)
        return linregress(self.x, logy)


class PolFit(Fit):
    """power law fit of form a*D**b"""
    def __init__(self, params=None, quess=(1, 0.2), **kwargs):
        super().__init__(params=params, name='polfit', **kwargs)
        self.quess = quess
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

    def find_fit(self, store_params=True, loglog=True, **kwargs):
        if not loglog:
            return super().find_fit(store_params=store_params, **kwargs)
        params, cov = self.find_fit_loglog()
        if store_params:
            self.params = params
            self.cov = cov
        return params, cov

    def find_fit_loglog(self):
        logx = np.log10(self.x)
        logy = np.log10(self.y)
        fitfunc = lambda p, x: p[0]+p[1]*x
        errfunc = lambda p, x, y: fitfunc(p, x)-y
        out = optimize.leastsq(errfunc, self.quess, args=(logx, logy),
                               full_output=True)
        pfinal = out[0]
        cov = out[1]
        a = 10**pfinal[0]
        b = pfinal[1]
        return (a, b), cov


gunn_kinzer = ExpFit(params=GUNN_KINZER)
set_plot_style()
