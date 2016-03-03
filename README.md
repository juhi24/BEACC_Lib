BAECC
=====

A Python 3 module for processing data from various (mostly in situ) precipitation instruments including Particle Imaging Package (PIP) and Pluvio^2 -pluviometers.

Installation
------------

The baecc-library requires python3 and compatible versions of numpy, scipy, pandas and matplotlib. 
The recommended way is to install these in a virtual enviroment.

You can create a virtual environment called py3env with virtualenv (python<=3.2)

    virtualenv -p /usr/bin/python3 py3env
    
OR pyvenv (python>=3.3)

    pyvenv py3env

Activate your virtual environment.

    source py3env/bin/activate
    
Then install the requirements inside the virtual environment using pip. Baecc requires `scipy`, `numpy`, `matplotlib`, `tstables` (`pytables`), `pytmatrix` and `pandas`.

    pip install PACKAGENAME

Clone the baecc python library.

    git clone https://github.com/juhi24/baecc.git
    
You should now have a everything ready to start working with the data.
  
Getting started
---------------

Files starting with `scr_` are example scripts.

You need PIP and pluviometer data either in raw ASCII or HDF5 format. 
Parts of the code can also utilize radar data using pytmatrix.

Data is selected using csv case configuration files. 
Minimum requirement for the configuration file is a header line with columns `start` and `end`, and at least one row (case) with corresponding datetime strings as shown here:

    start,end
    2014 21 February 21 UTC,2014 21 February 23 UTC

You can use different date formatting and specify a corresponding format string when initializing `EventsCollection` object.

Authors
-------

Jussi Tiira <jussi.tiira@helsinki.fi>

Davide Ori

Contributions also from Kalle Nordling

Full list of contributions can be viewed in the contributor tab in Github.

Licensing
---------

This work is licensed under GNU General Public v3 license. See `LICENSE` for more info.
