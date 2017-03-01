# coding: utf-8


class Cacher:
    """common methods to use msg cache"""
    def __init__(self, use_cache=True, storefilename='store.h5',
                 parent=None):
        self.use_cache = use_cache
        self.storefilename = storefilename
        self.parent = parent

    def msger(self, name, func, **kwargs):
        """Read from msgpack if caching is in use."""
        if self.use_cache:
            return self.msg_io(name, func, **kwargs)
        return func(**kwargs)

    def pickler(self, name, func, **kwargs):
        if self.use_cache:
            return self.pkl_io(name, func, **kwargs)
        return func(**kwargs)

    def cache_dir(self):
        """Return full path to cache directory."""
        if self.parent is None:
            return os.path.join(CACHE_DIR, self.fingerprint())
        return os.path.join(CACHE_DIR, self.parent.fingerprint(),
                            self.fingerprint())

    def store_path(self):
        """Return full path to hdf store file."""
        return os.path.join(self.cache_dir(), self.storefilename)

    def store_read(self, tablename, default_value=None, nocache_value=None):
        """Read from hdf store if using caching."""
        if self.use_cache:
            ensure_dir(self.cache_dir())
            try:
                with pd.HDFStore(self.store_path()) as store:
                    return store.get(tablename)
            except KeyError as err:
                warnings.warn("KeyError: {0} Using default value.".format(err))
                return default_value
        return nocache_value

    def store_write(self, tablename, data):
        """Write data to hdf store."""
        ensure_dir(self.cache_dir())
        with pd.HDFStore(self.store_path()) as store:
                store[tablename] = data

    def msg_io(self, name, func, **kwargs):
        """Read data from msgpack. If not available, calculate and store."""
        cd = self.cache_dir()
        msgpath = os.path.join(cd, name + MSGTLD)
        if os.path.isfile(msgpath):
            data = pd.read_msgpack(msgpath)
        else:
            ensure_dir(cd)
            data = func(**kwargs)
            data.to_msgpack(msgpath)
        return data

    def pkl_io(self, name, func, **kwargs):
        cd = self.cache_dir()
        pklpath = os.path.join(cd, name + PICKLETLD)
        if os.path.isfile(pklpath):
            with open(pklpath, 'rb') as cachefile:
                data = pickle.load(cachefile)
        else:
            ensure_dir(cd)
            data = func(**kwargs)
            with open(pklpath, 'wb') as cachefile:
                pickle.dump(data, cachefile, pickle.HIGHEST_PROTOCOL)
        return data

    def clear_cache(self, extra_files=None):
        """Remove cache files used by the Cacher object."""
        store = self.store_path()
        filelist = []
        if os.path.exists(store):
            filelist.append(store)
        if extra_files is not None:
            filelist.extend(extra_files)
        for f in filelist:
            os.remove(f)

    def fingerprint(self):
        """state-aware object identifier, immutable between sessions"""
        pass


class PrecipMeasurer:
    """parent for classes with precipitation measurements
    Either amount or acc (or both) methods should be overridden."""
    def __init__(self):
        pass

    def amount(self, **kwargs):
        """timestep precipitation in mm"""
        am = self.acc(**kwargs).diff()
        am.name = 'amount'
        return am

    def acc(self, **kwargs):
        """precipitation accumulation in mm"""
        acc = self.amount(**kwargs).cumsum()
        acc.name = 'accumulation'
        return acc

    def intensity(self, tdelta=None, **kwargs):
        """precipitation intensity in mm/h"""
        r = self.amount(**kwargs)
        if tdelta is None:
            tdelta = r.index.freq.delta
        frac = tdelta.apply(lambda t: np.timedelta64(1,'h')/t)
        intensity = frac * r
        intensity.name = 'intensity'
        return intensity


class InstrumentData(Cacher):
    """Parent for instrument data classes."""
    # TODO: Separate read_csv and __init__
    def __init__(self, filenames=None, data=None, hdf_table=None, use_cache=True):
        """Read from either ASCII data file or hdf5."""
        self.filenames = filenames
        if data is None:
            self.data = pd.DataFrame()
        else:
            self.data = data
        # if filtered data needed often, keep in memory
        self.stored_good_data = None    # set to None to disable
        if hdf_table is not None:
            self.name = hdf_table
            self.data = self.data.append(pd.read_hdf(filenames[0], hdf_table))
        Cacher.__init__(self, use_cache=use_cache)

    def __add__(self, other):
        combined = copy.deepcopy(self)
        combined.data = pd.concat([self.data, other.data])
        combined.clear_cache()
        return combined

    @classmethod
    def from_raw(cls, dtstrlist, subpath='', datadir=DATA_DIR):
        filelist = datafilelistloop(subpath, dtstrlist, datadir=datadir)
        return cls(filelist)

    def fingerprint(self):
        return fingerprint(str(self.data))

    def finish_init(self, dt_start, dt_end):
        """Sort and name index, cut time span."""
        self.data.sort_index(inplace=True)
        self.data.index.names = ['datetime']
        self.set_span(dt_start, dt_end)
        Cacher.__init__(self, storefilename=self.name + '.h5')

    def store_good_data(self, **kwargs):
        """Store good data to memory (to bypass recalculation of filters)."""
        self.stored_good_data = self.good_data(**kwargs)

    def parse_datetime(self):
        """Parse timestamps in data files. Used by class constructor."""
        pass

    def good_data(self):
        """Return useful data with filters and corrections applied."""
        if self.stored_good_data is not None:
            return self.stored_good_data
        return self.data

    def to_hdf(self, filename='../DATA/baecc.h5'):
        """Save object in hdf5 format."""
        self.data.to_hdf(filename, self.name, format='table', append=True)

    def between_datetime(self, date_start, date_end, inplace=False):
        """Limit the time span of data."""
        if inplace:
            instr = self
        else:
            instr = copy.deepcopy(self)
        instr.set_span(date_start, date_end)
        return instr

    def set_span(self, dt_start, dt_end):
        """Limit data to given time interval."""
        for dt in [dt_start, dt_end]:
            dt = pd.datetools.to_datetime(dt)
        self.data = self.data[dt_start:dt_end].copy()

    def grouped(self, varinterval=False, rule=None, col=None):
        if rule is None:
            rule = self.rule
        data = self.good_data()
        if col is not None:
            data = pd.DataFrame(data[col])
        if varinterval:
            grpd_data = pd.merge(data, rule, left_index=True, right_index=True)
            return grpd_data.groupby('group')
        return data.groupby(pd.Grouper(freq=rule, closed='right',
                                                   label='right'))