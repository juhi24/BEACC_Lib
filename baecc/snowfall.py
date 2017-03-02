"""Tools for estimating density and other properties of falling snow"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy
import warnings
import bisect
from itertools import cycle


def split_index(df, date=pd.datetime(2014,7,1), names=('first', 'second')):
    isfirst = df.index < date
    idf = pd.Series(isfirst, index=df.index)
    idf[isfirst] = names[0]
    idf[-isfirst] = names[1]
    tuples = list(zip(*(idf.values, idf.index.values)))
    index = pd.MultiIndex.from_tuples(tuples, names=('winter', 'datetime'))
    df.index = index
    return df


def before_after_col(df, date=pd.datetime(2014,7,1), colname='winter',
                     datecol=None):
    if datecol is None:
        dates = df.index
    else:
        dates = df[datecol]
    isfirst = dates > date
    df[colname] = isfirst.astype(int)
    return df


def remove_subplot_gaps(f, axis='row', axarr=None):
    if axarr is None:
            axarr=np.array(f.axes)
    adjust_kws = {}
    if axis=='row':
        if axarr.ndim==2:
            ax_whipe_ticks = axarr[:, 1:].flatten()
        else:
            ax_whipe_ticks = axarr[1:]
        adjust_kws['wspace'] = 0
        labels = [ax.get_yticklabels() for ax in ax_whipe_ticks]
    elif axis=='col':
        if axarr.ndim==2:
            ax_whipe_ticks = axarr[:-1, :].flatten()
        else:
            ax_whipe_ticks = axarr[:-1]
        adjust_kws['hspace'] = 0
        labels = [ax.get_xticklabels() for ax in ax_whipe_ticks]
    f.subplots_adjust(**adjust_kws)
    plt.setp(labels, visible=False)


def deprecation(message, stacklevel=2):
    """Issue DeprecationWarning"""
    warnings.warn(message, DeprecationWarning, stacklevel=stacklevel)


def daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days)):
        yield start_date + timedelta(n)


def scatterplot(x, y, c=None, kind='scatter', **kwargs):
    """scatter plot of two Series objects"""
    plotdata = read.merge_series(x, y)
    if c is not None:
        kwargs['c'] = c
    return plotdata.plot(kind=kind, x=x.name, y=y.name, **kwargs)


def combine_datasets(*datasets):
    eventslist = []
    for e in datasets:
        eventslist.append(e.events)
    combined = copy.deepcopy(datasets[0])
    combined.events = pd.concat(eventslist)
    combined.events.sort(columns='start', inplace=True)
    combined.events.reset_index(inplace=True)
    return combined


def plot_pairs(data, x='a', y='b', c=None, sizecol=None, scale=1,
               kind='scatter', groupby=None, ax=None, colorbar=False,
               markers='os^vD*p><', edgecolors='none', dtformat='%Y %b %d',
               split_date=None, **kwargs):
        """Easily plot parameters against each other."""
        if ax is None:
            ax = plt.gca()
        if c is not None:
            kwargs['c'] = c
        if sizecol is not None:
            kwargs['s'] = scale*np.sqrt(data[sizecol])
        if groupby is not None:
            groups = data.groupby(groupby)
            for (name, group), marker in zip(groups, cycle(markers)):
                colorbar = groups.case.first().iloc[0] == name and colorbar
                group.plot(ax=ax, x=x, y=y, marker=marker, kind=kind,
                           label=name, colorbar=colorbar,
                           edgecolors=edgecolors, **kwargs)
            return ax
        return data.plot(ax=ax, x=x, y=y, kind=kind, colorbar=colorbar,
                         edgecolors=edgecolors, **kwargs)


def find_interval(x, limits=(0, 100, 200, 1000)):
    """Find rightmost value less than x and leftmost value more than x."""
    i = bisect.bisect_right(limits, x)
    return limits[i-1:i+1]


def find_interval_df(s, limits):
    """Find intervals for Series s, output as a two-column DataFrame."""
    return s.apply(find_interval, limits=limits).apply(pd.Series)


def apply_rho_intervals(df, limits, rho_col='density'):
    """Add columns for density intervals."""
    data = df.copy()
    data[['rhomin', 'rhomax']] = find_interval_df(data[rho_col], limits)
    return data


def plot_vfit(case, dt, ax=None, extent=(0.3, 4, 0.5, 1.5),
              xtick_pos=(1, 2, 3, 4), tformat='%H:%M'):
    if ax is None:
        ax = plt.gca()
    vfit = case.instr['pipv'].fits.polfit[dt]
    vfit.plot(source_style='hex', unfiltered=True, ax=ax,
              source_kws={'gridsize': 20, 'extent': extent})
    ax.axis(extent)
    dt_start = case.instr['pluvio'].start_time()[dt]
    dt_end = pd.to_datetime(dt)
    ax.set_title('{0}–{1}'.format(dt_start.strftime(tformat),
                                  dt_end.strftime(tformat)))
    ax.set_xticks(xtick_pos)
    ax.tick_params(axis='both', direction='out', length=4)
    ax.legend()


def series_cdf(series):
    """CDF from a pandas Series"""
    srtd = series.sort_values()
    cum_dist = np.linspace(0, 1, series.size)
    return pd.Series(cum_dist, index=srtd)
