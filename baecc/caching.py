# coding: utf-8
import os
import warnings
import pickle
import hashlib
import pandas as pd
from j24 import home, ensure_dir

MSGTLD = '.msg'
PICKLETLD = '.pkl'
CACHE_DIR = os.path.join(home(), '.cache', 'baecc')


def fingerprint(string):
    return hashlib.sha256(string.encode('utf-8')).hexdigest()[-12:]


def hash_dict(d):
    return fingerprint(str(sorted(d.items())))


def combine2str(*identifiers):
    return ''.join(tuple(map(str, identifiers)))


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