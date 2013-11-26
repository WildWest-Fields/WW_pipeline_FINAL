#! /usr/bin/env python

'''
ABOUT:
This script runs tweakback to apply the WCS solution from a drizzled image (since aligned to a reference)
to its constituent flt images as determined from the header of the drizzled image. 
The flt images may then be drizzled to any footprint.

DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER, 2013

HISTORY:
September 2013: Original script (v0.1).

FUTURE IMPROVEMENTS:

USE:
python run_tweakback.py	-drz drz.fits	--> apply WCS solution from drz.fits to constituent flts (must be located in cwd).

'''

__author__='D.M. HAMMER'
__version__= 0.1

import os, glob, argparse, pyfits, pdb
from pyraf import iraf
import numpy as np
from drizzlepac import tweakback
from stsci.tools import teal


if __name__=='__main__':
    # -- Parse input parameters
    parser = argparse.ArgumentParser(description='Run tweakback to apply WCS from drz to fl? images.')
    parser.add_argument('-drz', '--drzim',default='NONE', type=str, help='Input drizzled image (no wildcards). \
    				Default is NONE - requires input.')
    options = parser.parse_args()

    # -- initialize names of drizzled/flt(c) images
    if options.drzim == 'NONE': raise Exception('Must input a drizzled image.')
    drzim = options.drzim

    # -- run tweakback
    iraf.unlearn('tweakback')
    teal.unlearn('tweakback')
    tweakback.tweakback(drzim,verbose=True)
