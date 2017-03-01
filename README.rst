baecc
=====

A Python 3 module for processing data from various (mostly in situ) precipitation instruments including Particle Imaging Package (PIP) and Pluvio^2 -pluviometers.

Requirements
------------

baecc requires ``numpy``, ``scipy``, ``pandas``, ``tstables``, ``pytmatrix``, ``netCDF4`` and ``seaborn``

Installation
------------

It is recommended to install the module in a virtual python3 environment (e.g. using anaconda). To install the module, enter your virtual environment and run

::

    python setup.py install

Getting started
---------------

There are some example scripts in the scripts directory.

You need PIP and pluviometer data either in raw ASCII or HDF5 (recommended) format.
Parts of the code can also utilize radar data with pytmatrix.
You can convert your data to a HDF5 database using ``read.batch_create_hdf``.

Data is selected using csv case configuration files.
Minimum requirement for the configuration file is a header line with columns ``start`` and ``end``, and at least one row (case) with corresponding datetime strings as shown here:

::

    start,end
    2014 21 February 21 UTC,2014 21 February 23 UTC

You can use different date formatting and specify a corresponding format string when initializing ``EventsCollection`` object.

Authors
-------

Jussi Tiira <jussi.tiira@helsinki.fi>

Davide Ori

Contributions also from Kalle Nordling

Full list of contributions can be viewed in the contributor tab in Github.

Licensing
---------

This work is licensed under GNU General Public v3 license. See ``LICENSE`` for more info.
