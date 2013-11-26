#! /usr/bin/env python

'''
ABOUT:
This program runs tweakreg to align drizzled images to a reference image using external catalogs.

DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER, 2013

HISTORY:
July 2013: Original script (v0.1).

FUTURE IMPROVEMENTS:

USE:
python run_tweakreg_drz.py
'''

__author__='D.M. HAMMER'
__version__= 0.1

import os, glob, argparse, pyfits, pdb
import numpy as np
from drizzlepac import tweakreg
from stsci.tools import teal


if __name__=='__main__':
    # Parse input parameters
    parser = argparse.ArgumentParser(description='Run tweakreg on input images using custom source catalogs.')
    parser.add_argument('-im',  '--images',default='NONE', type=str, help='Input drizzled image(s). \
				 Default is NONE - must input images.')
    parser.add_argument('-cf',  '--catfile',default='catfile.sex', type=str, help='Input name of catfile. \
				 Default is catfile.sex.')
    parser.add_argument('-log', '--logfile',default='tweakreg.log', type=str, help='Input name of tweakreg log file. \
				 Default is tweakreg.log')
    parser.add_argument('-rim', '--refim',default='', type=str, help='Input reference image. \
                                 There is no default - must be entered.')
    parser.add_argument('-rcat','--refcat',default='', type=str, help='Input image file(s). \
                                 There is no default - must be entered.')
    parser.add_argument('-xc',  '--xcol',default=8, type=int, help='Input catalog column that corresponds to x. \
                                 Default = SExtractor XWIN column.')
    parser.add_argument('-yc',  '--ycol',default=9, type=int, help='Input catalog column that corresponds to y. \
                                 Default = SExtractor YWIN column.')
    parser.add_argument('-rxc', '--refxcol',default=8, type=int, help='Input (ref)catalog column that corresponds to x. \
                                 Default = SExtractor XWIN column.')
    parser.add_argument('-ryc', '--refycol',default=9, type=int, help='Input (ref)catalog column that corresponds to y. \
                                 Default = SExtractor YWIN column.')
    parser.add_argument('-runit','--refunit',default='pixels', type=str, help='Input tweakreg "refxyunits" parameter. \
                                  Default = "pixels".')
    parser.add_argument('-fit', '--fittype',default='rscale', type=str, help='Input tweakreg "fitgeometry" parameter. \
                                 Default = "rscale", i.e., shift, rotation, and pixel size.')
    				 
    options = parser.parse_args()
    logfile = options.logfile
    catfilename = options.catfile
    xcol = options.xcol
    ycol = options.ycol
    rxcol = options.refxcol
    rycol = options.refycol
    runit = options.refunit
    fittype = options.fittype

    # -- initialize image list and print to file
    if options.images == 'NONE': raise Exception('Must input drizzled images.')
    im = glob.glob(options.images)
    im.sort()
    f = open('imlist.dat', 'w')
    for ff in im: f.write(ff+'\n')
    f.close()
    
    # -- initialize reference image/catalog info
    irefim = glob.glob(options.refim)
    irefim.sort()
    if ((options.refim != '') & (len(irefim) == 0)):
        raise Exception('Could not find reference images matching '+options.refim+' in cwd.') 
    irefcat = glob.glob(options.refcat)
    irefcat.sort()
    if ((options.refcat != '') & (len(irefcat) == 0)):
        raise Exception('Could not find reference catalog matching '+options.refcat+' in cwd.')

    if len(irefim) > 0:
        irefim = irefim[0]
        irefcat = irefcat[0]
	USE_REF = True
    else:
        irefim = ''
        irefcat = ''
        USE_REF = False


    # -- set conv width based on filter name (IR = 2.5; optical=3.5)
    fheader = pyfits.getheader(im[0])
    instr = fheader['INSTRUME']
    if instr == 'WFC3':
        filtname = fheader['FILTER']
    elif instr == 'ACS':
        filtname = fheader['FILTER1']
        if filtname[0] == 'C': filtname = fheader['FILTER2']  
    else: raise Exception('Instrument '+instr+' not covered in our case list.')

    if filtname[0:2] == 'F1': conv_wid = 2.5
    else: conv_wid = 3.5


    # -- run tweakreg (SExtractor catalog: XWIN/YWIN cols 8/9; X_IMAGE/Y_IMAGE cols 2/3)
    teal.unlearn('tweakreg')

    if USE_REF:
        tweakreg.TweakReg('@imlist.dat',runfile=logfile,shiftfile=True,outshifts='shift.dat',updatehdr=True, \
			fitgeometry=fittype,catfile=catfilename,xcol=xcol,ycol=ycol, \
			refimage=irefim,refcat=irefcat,refxyunits=runit,refxcol=rxcol,refycol=rycol, \
			conv_width=conv_wid,searchrad=3.0,nclip=7,see2dplot=False,residplot='No plot')
                        #headerlet=True,attach=False,hdrname='TWEAKDRZ')
    else:
        tweakreg.TweakReg('@imlist.dat',runfile=logfile,shiftfile=True,outshifts='shift.dat',updatehdr=True, \
			fitgeometry=fittype,catfile=catfilename,xcol=xcol,ycol=ycol, \
			conv_width=conv_wid,searchrad=3.0,nclip=7,see2dplot=False,residplot='No plot')


    # -- remove unwanted tweakreg files
    tmp=np.concatenate((glob.glob('*.coo'),glob.glob('*catalog.match'),glob.glob('shift*.fits')))
    for ff in tmp: os.remove(ff)

