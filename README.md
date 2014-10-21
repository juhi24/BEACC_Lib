BAECC
=====

A Python 3 module to work with baecc campaign in situ data

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
    
Then install the requirements inside the virtual environment using pip. Baecc requires `numpy`, `matplotlib` and `pandas`.

    pip install PACKAGENAME

Clone the baecc python library.

    git clone https://github.com/juhi24/baecc.git
    
You should now have a everything ready to start working with the data.
  
Getting started
---------------

Files starting with `scr_` are example scripts.
