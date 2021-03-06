# coding: utf-8
from __future__ import absolute_import, division, print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
from os import path
from scipy import stats
from scipy.optimize import fmin, minimize
import baecc
from baecc import fit, instruments, caching
from j24 import ensure_dir

SUBPATH = 'PIP/a_Velocity_Tables/004%s/*2.dat'

def kde(x, y):
    values = np.vstack((x, y))
    return stats.gaussian_kde(values)


def bindata(xmin, xmax, data, ymin=None, ymax=None, xname='x', yname='y'):
    """Return data that falls into given x bin."""
    cond = '%s > %s and %s < %s' % (xname, xmin, xname, xmax)
    if ymin is not None and ymax is not None:
        cond += ' and %s > %s and %s < %s' % (yname, ymin, yname, ymax)
    return data.query(cond)


def filter_outlier(X, Y, Z, data, xname='x', yname='y', frac=0.5,
                   bin_limit_multiplier=0.05):
    filtered = pd.DataFrame()
    std = []
    x = X[0, :]
    y = Y[:, 0]
    xwidth = (x[-1]-x[0])/len(x)    # TODO: check if correct
    #print('xw: %s' % xwidth)
    xlims = np.append(x-0.5*xwidth, x[-1]+0.5*xwidth)
    xmin = xlims[:-1]
    xmax = xlims[1:]
    ymin = []
    ymax = []
    data_count = data.count()[0]
    for i in range(0, Z.shape[1]):  # loop through bins
        z = Z[:, i]
        z_lim = z.max()*frac
        y_fltr = y[z > z_lim]       # FWHM when frac=0.5
        if y_fltr.size == 0:
            ymin.append(None)
            ymax.append(None)
            continue
        ymin.append(y_fltr[0])
        ymax.append(y_fltr[-1])
        bin_data = bindata(xmin[i], xmax[i], data, ymin=ymin[-1], ymax=ymax[-1],
                           xname=xname, yname=yname)
        bin_count = bin_data.count()[0]
        count_limit = bin_limit_multiplier*xwidth*data_count
        if bin_count < count_limit:
            ymin[-1] = None
            ymax[-1] = None
            continue
        std.append(bin_data[yname].std())
        filtered = filtered.append(bin_data)
    # return filtered data, stds and half width at frac*kde_max
    return filtered, np.array(std), xlims, ymin, ymax


class PipV(instruments.InstrumentData):
    """PIP particle velocity and diameter data handling"""
    def __init__(self, filenames=None, dt_start=None, dt_end=None, **kwargs):
        """Create a PipV object using data from a list of PIP velocity table
        files."""
        instruments.InstrumentData.__init__(self, filenames, **kwargs)
        self.name = 'pip_vel'
        self.dmin = 0.3   # shortest Wad_Dia where data is good
        self.dmax = 25.8
        self.vmin = 0.5
        self.d_col = 'd_voleq' # equivalent volume
        #self.d_col = 'Wad_Dia' # equivalent area
        self._fits = pd.DataFrame()
        # num=511 --> binwidth 0.05
        if baecc.DEBUG:
            num = 103
        else:
            num = 409
        self.dbins = np.linspace(self.dmin, self.dmax, num=num)
        self._std = pd.DataFrame(columns=self.dbins)
        # half width at fraction of maximum
        self._hwfm = pd.DataFrame(columns=self.dbins)
        self.default_fit = fit.PolFit
        self.flip = False
        self.loglog = True # use loglog method with power law fitting
        if self.data.empty:
            print('Reading PIP particle velocity data...')
            for filename in filenames:
                print('.', end='')
                self.current_file = filename
                if int(filename[-23:-15]) > 20141124:
                    newdata = pd.read_csv(filename,
                                          engine='python', sep='\t',
                                          skipinitialspace=True, skiprows=8,
                                          skip_footer=1,
                                          parse_dates={'datetime':['minute_p']},
                                          date_parser=self.parse_datetime,
                                          verbose=baecc.DEBUG)
                else:
                    newdata = pd.read_csv(filename,
                                          delim_whitespace=True,
                                          skiprows=8,
                                          parse_dates={'datetime':['minute_p']},
                                          date_parser=self.parse_datetime,
                                          verbose=baecc.DEBUG)
                newdata = newdata[newdata['RecNum']>-99]
                if not newdata.empty:
                    newdata.rename_axis({'vel_v_1': 'vel_v',
                                         'vel_h_1': 'vel_h',
                                         'vel_v_2': 'vel_v',
                                         'vel_h_2': 'vel_h'},
                                        axis=1, inplace=True)
                    newdata.set_index(['datetime', 'Part_ID', 'RecNum'],
                                      inplace=True)
                    g = newdata.groupby(level=['datetime', 'Part_ID'])
                    newdata = g.mean()
                    self.data = self.data.append(newdata)
            print()
            if len(self.data.index):
                self.data = self.data[self.data.vel_v.notnull()]
            self.data.reset_index(level=1, inplace=True)
            self.data = self.data.astype(float)
        self.finish_init(dt_start, dt_end)

    @property
    def rule(self):
        if self.fits.empty:
            return None
        #return self.fits.index.freq.freqstr
        return self.fits.index.freqstr

    @property
    def binwidth(self):
        d = self.dbins
        return (d[-1]-d[0])/(len(d)-1)

    @property
    def fits(self):
        return self.store_read('fits', default_value=pd.DataFrame(),
                               nocache_value=self._fits)

    @fits.setter
    def fits(self, fits):
        if self.use_cache:
            self.store_write('fits', fits)
        else:
            self._fits = fits

    @property
    def std(self):
        return self.store_read('std',
                               default_value=pd.DataFrame(columns=self.dbins),
                               nocache_value=self._std)

    @std.setter
    def std(self, std):
        if self.use_cache:
            self.store_write('std', std)
        else:
            self._std = std

    @property
    def hwfm(self):
        return self.store_read('hwfm',
                               default_value=pd.DataFrame(columns=self.dbins),
                               nocache_value=self._hwfm)

    @hwfm.setter
    def hwfm(self, hwfm):
        if self.use_cache:
            self.store_write('hwfm', hwfm)
        else:
            self._hwfm = hwfm

    @classmethod
    def from_raw(cls, *args, subpath=SUBPATH, **kwargs):
        return super().from_raw(*args, subpath=subpath, **kwargs)

    def fingerprint(self):
        identifiers = (super().fingerprint(), self.flip, self.dbins, self.loglog)
        idstr = caching.combine2str(*identifiers)
        return caching.fingerprint(idstr)

    def v(self, d, fitclass=None, varinterval=True, rule=None):
        """velocities according to fits for given diameter"""
        if fitclass is None:
            fitclass = self.default_fit
        if rule is None:
            rule = self.rule
        if self.fits.empty:
            self.find_fits(rule, fitclass=fitclass, varinterval=varinterval)
        elif not varinterval:
            if pd.datetools.to_offset(rule) != self.fits.index.freq:
                print('different sampling freq')
                self.find_fits(rule, fitclass=fitclass,
                               varinterval=varinterval)
        v = []
        for vfit in self.fits[fitclass.name].values:
            v.append(vfit.func(d))
        return pd.Series(v, index=self.fits.index, name='v')

    def lwc(self, rule='1min'):
        """liquid water content"""
        d3 = self.good_data()[self.d_col]**3
        return d3.resample(rule, how=np.sum, closed='right', label='right')

    def parse_datetime(self, mm):
        datestr = self.current_file.split('/')[-1].split('_')[0]
        yr = int(datestr[3:7])
        mo = int(datestr[7:9])
        dd = int(datestr[9:11])
        hh = int(datestr[11:13])
        return datetime.datetime(yr, mo, dd, hh, int(mm))#,
                                 #tzinfo=datetime.timezone.utc)

    def good_data(self):
        if self.stored_good_data is not None:
            return self.stored_good_data
        query_str = 'Wad_Dia > {0} & vel_v > {1}'.format(self.dmin, self.vmin)
        data = self.data.query(query_str)
        data.loc[:,'d_voleq'] = data.Wad_Dia/instruments.PHI
        return data

    def filter_outlier(self, data, frac=0.5, flip=False):
        """Filter outliers using KDE"""
        if data is None:
            data = self.good_data()
        if flip:
            X, Y, Z = self.kde_grid(data)
            return filter_outlier(Y.T, X.T, Z.T, data, xname='vel_v',
                                  frac=frac, yname=self.d_col)
        return filter_outlier(*self.kde_grid(data), data=data,
                              xname=self.d_col, yname='vel_v', frac=frac)

    def frac_larger(self, d):
        """Return fraction of particles that are larger than d."""
        vdata = self.good_data()
        d_data = vdata[self.d_col]
        return vdata[d_data > d].vel_v.count()/vdata[d_data < d].vel_v.count()

    def d_cut(self, frac=0.05, d_guess=2):
        """Return d for which given fraction of particles are larger."""
        dcost = lambda d: abs(self.frac_larger(d[0])-frac)
        return fmin(dcost, d_guess)[0]

    def find_fit(self, fitclass=None, data=None, use_kde_peak=False,
                 cut_d=False, frac=0.5, use_curve_fit=True, bin_num_min=5,
                 filter_outliers=True, name=None, try_flip=None,
                 plot_flip=False, force_flip=False, cut_kws={}, **kwargs):
        """Find and store a fit for either raw data or kde."""
        # TODO: clean this mess
        def too_few_particles(use_curve_fit, use_kde_peak):
            print('Too few particles.')
            return False, False
        std = pd.DataFrame()
        hwfm = pd.DataFrame()
        if try_flip is None:
            try_flip = self.flip
        if force_flip:
            try_flip = True
        if data is None:
            data = self.good_data()
        origdata = data.copy() # do not rewrite
        if fitclass is None:
            fitclass = self.default_fit
        vfit = fitclass()
        partcount = data.count()[0]
        # increased from 5 to 10
        if partcount < 10 and (use_curve_fit or use_kde_peak):
            use_curve_fit, use_kde_peak = too_few_particles(use_curve_fit,
                                                            use_kde_peak)
        elif filter_outliers:
            data, stdarr, xlims, ymin, ymax = self.filter_outlier(data=data,
                                                             frac=frac)
            vfit.fltr_upper_x = xlims
            vfit.fltr_lower_x = xlims
            vfit.fltr_upper_y = ymax + [ymax[-1]]
            vfit.fltr_lower_y = ymin + [ymin[-1]]
            # TODO: Make a proper conditional to include this when needen, if needed
            if False:
                datao = self.filter_outlier(data=data, frac=frac, flip=True)[0]
                fltrcount = datao.count()[0]
                if fltrcount < 2 and (use_curve_fit or use_kde_peak):
                    use_curve_fit, use_kde_peak = too_few_particles(use_curve_fit,
                                                                    use_kde_peak)
            if name is not None:
                stdarr.resize(self.dbins.size, refcheck=False)
                std = pd.DataFrame(pd.Series(stdarr,
                                             index=self.dbins, name=name)).T
                # TODO: empty hwfm
                hwfm = pd.DataFrame(pd.Series(index=self.dbins, name=name)).T
                for df in [std, hwfm]:
                    df.index.name = 'datetime'
        else:
            print('Did not apply filter.')
        if use_kde_peak:
            d, v = self.kde_peak(data=data)
        else:
            d = data[self.d_col].values
            v = data.vel_v.values
        if cut_d:
            dcut = self.d_cut(**cut_kws)
            d = d[d < dcut]
            v = v[d < dcut]
        if use_kde_peak:
            num = np.array([bindata(diam, self.binwidth, data=data,
                                    xname=self.d_col,
                                    yname='vel_v').vel_v.count() for diam in d])
            d = d[num > bin_num_min]
            v = v[num > bin_num_min]
            sig = [self.dbin(diam, self.binwidth, data=data, xname=self.d_col,
                             yname='vel_v').vel_v.sem() for diam in d]
        else:
            sig = np.ones(d.size)
        vfit.x = d
        vfit.y = v
        vfit.x_unfiltered = origdata[self.d_col].values
        vfit.y_unfiltered = origdata.vel_v.values
        if use_curve_fit:
            unflipped_kws = kwargs
            if self.default_fit.name == fit.PolFit.name:
                unflipped_kws['loglog'] = self.loglog
            params, pcov = vfit.find_fit(**unflipped_kws)
            perr = vfit.perr()   # standard errors of d, v
            if try_flip and not use_kde_peak:
                fiti = fitclass(flipped=True)
                datai, stdarri, xlimsi, ymini, ymaxi = self.filter_outlier(data=data,
                                                               frac=frac,
                                                               flip=True)
                fiti.fltr_upper_y = xlimsi
                fiti.fltr_lower_y = xlimsi
                fiti.fltr_upper_x = ymaxi + [ymaxi[-1]]
                fiti.fltr_lower_x = ymini + [ymini[-1]]
                fiti.x = datai[self.d_col].values
                fiti.y = datai.vel_v.values
                fiti.x_unfiltered = origdata[self.d_col].values
                fiti.y_unfiltered = origdata.vel_v.values
                paramsi, pcovi = fiti.find_fit(**kwargs)
                perri = fiti.perr()
                if plot_flip:
                    f, axarr = plt.subplots(1, 3, sharey=True, sharex=True,
                                            figsize=(12, 6))
                    for ax in axarr:
                        fiti.plot(ax=ax, label='flipped %.4f' % perri[1])
        else:
            result = minimize(vfit.cost, vfit.quess, method='Nelder-Mead',
                              args=(d, v, sig))
            vfit.params = result.x
        if plot_flip and use_curve_fit:
            self.plot_flip(axarr, f, vfit, data, datai, origdata, perr)
        desfmt = '{0:.4f}'
        fitstr = 'standard'
        errstr = ''
        fitout = vfit
        if use_curve_fit:
            errstr += 'err: std ' + desfmt.format(perr[1])
            if try_flip:
                errstr += ', inv ' + desfmt.format(perri[1])
                if (perr[1] > perri[1] and fiti.is_good()) or force_flip:
                    fitout = fiti
                    fitstr = 'flipped'
        print(fitstr + ' fit: ' + str(fitout) + '; ' + str(partcount) + ' particles; ' + errstr)
        return fitout, std, hwfm    # TODO: wrong std, hwfm when flipped

    def plot_flip(self, axarr, f, vfit, data, datai, origdata, perr):
        filterstr = ['D-binned filter', 'v-binned filter', 'unfiltered']
        for ax in axarr:
            vfit.plot(ax=ax, label='original %.4f' % perr[1])
        self.plot(ax=axarr[0], data=data, ymax=3)
        self.plot(ax=axarr[1], data=datai, ymax=3)
        self.plot(ax=axarr[2], data=origdata, ymax=3)
        for i, ax in enumerate(axarr):
            ax.set_title(ax.get_title() + ', ' + filterstr[i])
        plt.legend()
        f.tight_layout()
        fname = data.index[-1].strftime('%H%M.eps')
        datedir = data.index[-1].strftime('%Y%m%d')
        f.savefig(path.join(ensure_dir(path.join(baecc.RESULTS_DIR, 'pip2015',
                                                       'fitcomparison',
                                                       datedir)), fname))
        return axarr

    def find_fits(self, rule, varinterval=True, empty_on_fail=True, **kwargs):
        print('Calculating velocity fits for given sampling frequency...')
        names = []
        fits = []
        stds = []
        hwfms = []
        for name, group in self.grouped(rule=rule, varinterval=varinterval):
            try:
                newfit, std, hwfm = self.find_fit(data=group, name=name,
                                                  try_flip=self.flip,
                                                  **kwargs)
            except RuntimeError as err:
                print('%s: %s' % (name, err))
                print('Particle count: %s' % group.vel_v.count())
                if len(fits) == 0 or empty_on_fail:
                    print('Using an empty fit')
                    newfit = self.default_fit()
                    std = np.nan
                    hwfm = np.nan
                else:
                    print('Using fit from previous time step.')
                    newfit = fits[-1]
                    std = stds[-1]
                    hwfm = hwfms[-1]
            fits.append(newfit)
            names.append(name)
            stds.append(std)
            hwfms.append(hwfm)
        self.std = pd.concat(stds)
        self.hwfm = pd.concat(hwfms)
        if varinterval:
            timestamps = names
        else:
            timestamps = pd.DatetimeIndex(names, freq=rule)
        if self.fits.empty:
            self.fits = pd.DataFrame(fits, index=timestamps,
                                     columns=[newfit.name])
        elif self.fits.index.equals(timestamps):
            self.fits[newfit.name] = fits
        else:
            self.fits = pd.DataFrame(fits, index=timestamps,
                                     columns=[newfit.name])
        return self.fits

    def fit_params(self, fit_type=None):
        """Return DataFrame of fit parameters."""
        if fit_type is None:
            fit_type = self.default_fit.name
        letter = 'abcdef'
        params = self.fits[fit_type].apply(lambda vfit: vfit.params)
        paramlist = []
        for i in range(len(params.values[0])):
            param = params.apply(lambda p: p[i])
            param.name = letter[i]
            paramlist.append(param)
        return baecc.merge_multiseries(*paramlist) # TODO: replace with concat

    def partcount(self, rule, varinterval):
        return self.grouped(rule=rule, varinterval=varinterval).Part_ID.count()

    def kde(self, data=None):
        """kernel-density estimate of d,v data using gaussian kernels"""
        if data is None:
            data = self.good_data()
        d = data[self.d_col].values
        v = data.vel_v.values
        return kde(d, v)

    def grids(self, data=None):
        if data is None:
            data = self.good_data()
        d = data[self.d_col].values
        v = data.vel_v.values
        dmax = d.max()+20*self.binwidth
        dbins = self.dbins[self.dbins < dmax]
        num_vbins = round(len(self.dbins)/5)
        return np.meshgrid(dbins, np.linspace(v.min(), v.max(), num_vbins))

    def kde_grid(self, data=None):
        """Calculate kernel-density estimate with given resolution."""
        X, Y = self.grids(data)
        points = np.vstack([X.ravel(), Y.ravel()])
        kernel = self.kde(data)
        Z = np.reshape(kernel(points).T, X.shape)
        return X, Y, Z

    def kde_peak(self, **kwargs):
        """the most propable velocities for each diameter in grid"""
        D, V, Z = self.kde_grid(**kwargs)
        x = D[0, :]
        y = V[:, 0][Z.argmax(axis=0)]
        return x, y

    def plot_kde(self, ax=None, **kwargs):
        """Plot kde grid."""
        if ax is None:
            ax = plt.gca()
        D, V, Z = self.kde_grid(**kwargs)
        pc = ax.pcolor(D, V, Z, cmap=plt.cm.gist_earth_r)
        return pc

    def plot_fits(self, fit_type=None, savefigs=False, path='.',
                  fname='%Y%m%d_%H%M.eps', **kwargs):
        if fit_type is None:
            fit_type = self.default_fit.name    # as a string
        fits = self.fits[fit_type]
        axlist = []
        flist = []
        for i, vfit in fits.iteritems():
            flist.append(plt.figure())
            axlist.append(vfit.plot(**kwargs))
            plt.savefig(path.join(path, i.strftime(fname)))
        return pd.DataFrame(index=fits.index, data={'fit':fits, 'fig':flist,
                                                    'ax':axlist})

    def plot_fit(self, tstep=None, **kwargs):
        if tstep is None:
            fits = [self.find_fit()[0]]
        else:
            fits = self.fits.loc[tstep].values
        for vfit in fits:
            vfit.plot(**kwargs)

    def plot(self, data=None, hexbin=True, ax=None, xmax=None, ymax=None,
             show_particle_count=False, colormap='gray_r', ygrid=True,
             hexsize=12, **kwargs):
        """Plot velocity data."""
        if ax is None:
            ax = plt.gca()
        if data is None:
            data = self.good_data()
        margin = 0.1
        if xmax is None:
            xmax = np.ceil(self.d_cut(frac=0.05))
        right = xmax-2+margin
        partcount = data.Part_ID.count()
        if partcount < 1:
            return ax
        if hexbin:
            data.plot(x=self.d_col, y='vel_v', kind='hexbin', label='hexbinned',
                      ax=ax, gridsize=int(hexsize*data[self.d_col].max()**0.5),
                      colormap=colormap, **kwargs)
        else:
            data.plot(x=self.d_col, y='vel_v', style=',', ax=ax,
                      alpha=0.2, color='black', label='pip raw', **kwargs)
        #fit.gunn_kinzer.plot(dmax=20, label='Gunn&Kinzer', ax=ax, zorder=5, ls='--')
        if ymax is None:
            ymax = data.vel_v.max() + margin
        ax.axis([0, xmax, 0, ymax])
        ax.yaxis.grid(ygrid)
        t_start = data.index[0]-datetime.timedelta(minutes=1)
        t_end = data.index[-1]
        label_format = '%H:%M'
        ax.set_title('%s-%s UTC' % (t_start.strftime(label_format),
                                    t_end.strftime(label_format)))
        if show_particle_count:
            ax.text(right, margin, 'particle count: %s' % str(partcount))
        ax.set_ylabel('Fall velocity (m/s)')
        ax.set_xlabel('D (mm)')
        ax.legend(loc='upper right')
        return ax

    def plots(self, rule=None, separate=True, peak=False, save=False, ncols=1,
              prefix='', suffix='.png', ymax=None, plotfit=True, savedir=None,
              **kwargs):
        """Plot datapoints and fit for each timestep."""
        ngroups = self.grouped(rule=rule).ngroups
        #nrows = int(np.ceil(ngroups/ncols))
        if not separate:
            f, axarr = plt.subplots(1, ngroups, sharex='col', sharey='row',
                                    figsize=(ngroups*8, 7), tight_layout=True)
        else:
            axarr = []
            farr = []
        if ymax is None:
            self.good_data().vel_v.max()
        for i, (name, group) in enumerate(self.grouped(rule=rule)):
            if separate:
                f = plt.figure(dpi=175, figsize=(3.5, 3))
                ax = plt.gca()
                farr.append(f)
                axarr.append(ax)
            if group.Part_ID.count() < 1:
                continue
            if plotfit:
                self.plot_fit(tstep=name, zorder=6, ax=axarr[i], marker=',',
                              alpha=0.3)
            self.plot(data=group, ax=axarr[i],
                      ymax=ymax, **kwargs)
            f.tight_layout()
            if peak:
                axarr[i].scatter(*self.kde_peak(data=group), label='kde peak')
        return axarr
