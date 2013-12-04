#! /usr/bin/env python

'''
ABOUT:
This program drizzles a stack of flt images.

DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER, 2013

HISTORY:
July 2013: Original script (v0.1).
Oct. 2013: Removed iraf "unlearn" (unnecessary).


FUTURE IMPROVEMENTS:

USE:
python run_drizzle.py -fk square                -->drizzle all *fl?.fits images in cwd using final_kernel='square'.
python run_drizzle.py -im i*fl?.fits -c 12      -->drizzle only "i*fl?.fits" files in cwd using 12 cores.
python run_drizzle.py -pf 0.5 -ps 0.04          -->same as above but use final_pixfrac=0.5 & final_scale=0.04".
python run_drizzle.py -fb 96 -wht EXP           -->same as above but use final_bits=32,64 & final_wht_type='EXP'.
python run_drizzle.py -ra 196. -dec 27.0 -nx 8000 -ny 8000 -->"" but drizzle to footprint w/ra,dec center & 8000x8000 pix.

'''

__author__='D.M. HAMMER'
__version__= 0.1

import os, glob, argparse, pyfits, pdb
import numpy as np
from drizzlepac import astrodrizzle
from stsci.tools import teal


if __name__=='__main__':
    # Parse input parameters
    parser = argparse.ArgumentParser(description='Run tweakreg on input images using custom source catalogs.')
    parser.add_argument('-im',  '--images',default='*fl?.fits', type=str, help='Input image file(s). \
                                 Default is all _fl? images in current directory.')
    parser.add_argument('-out', '--outname',default='NA', type=str,help='Prefix for AD output.')
    parser.add_argument('-pf',  '--pixfrac',default=0.6, type=float, help='Input AstroDrizzle "final_pixfrac" parameter. \
                                 Default is 0.6 - expected value for large ACS stacks.')
    parser.add_argument('-ps',  '--pixscale',default=0.03, type=float, help='Input AstroDrizzle "final_scale" parameter. \
                                 Default value is 0.03" - expected value for large ACS stacks.')
    parser.add_argument('-fb',  '--finalbits',default=0, type=int, help='Input AstroDrizzle "final_bits" parameter. \
                                 Default value is 0 (assumes several visits). \
                                 A single visit drizzle should adopt 64 or 96 for WFC3 IR and ACS, respectively.')
    parser.add_argument('-fk',  '--finalkernel',default='gaussian', type=str, choices=['square','point','gaussian',\
                                 'turbo','tophat','lanczos3'],help='Input AstroDrizzle "final_kernel" parameter. \
                                 Default value is "gaussian".')
    parser.add_argument('-wht', '--whttype',default='IVM', type=str, choices=['EXP','ERR','IVM'],help='Input AstroDrizzle \
                                 "final_wht_type" parameter. Default value is IVM (source-free inv. variance map.')
    parser.add_argument('-ra',  '--finalra',default=None, type=float, help='Input AstroDrizzle "final_ra" parameter. \
                                 Default value indicates no footprint (AD selects image dimensions).')
    parser.add_argument('-dec', '--finaldec',default=None, type=float, help='Input AstroDrizzle "final_dec" parameter. \
                                 Default value indicates no footprint (AD selects image dimensions).')
    parser.add_argument('-nx',  '--finaloutnx',default=None, type=int, help='Input AstroDrizzle "final_outnx" parameter. \
                                 Default value indicates no restrictions on no. columns in output image.')
    parser.add_argument('-ny',  '--finaloutny',default=None, type=int, help='Input AstroDrizzle "final_outny" parameter. \
                                 Default value indicates no restrictions on no. rows in output image.')
    parser.add_argument('-c',   '--cores',default=4, type=int, help='Input number of processing cores used by AD. \
                                 Default is 4.')
    options = parser.parse_args()


    # -- initialize variables to input arguments
    NCORES = options.cores
    outname = options.outname
    pixscale = options.pixscale
    pixfrac = options.pixfrac
    finalbits = options.finalbits
    finalwht = options.whttype.upper()
    finalkernel = options.finalkernel.lower()
    finalra = options.finalra
    finaldec = options.finaldec
    finalnx = options.finaloutnx
    finalny = options.finaloutny


    # -- print image list to file ("imlist.dat")
    imlist = glob.glob(options.images)
    imlist.sort()
    f = open('imlist.dat', 'w')
    for ff in imlist: f.write(ff+'\n')
    f.close()

    # -- get instrument and filter name
    instrum = pyfits.getheader(imlist[0])['INSTRUME']
    if instrum == 'WFC3': filtname = pyfits.getheader(imlist[0])['FILTER']
    elif instrum == 'ACS':
        filtname = pyfits.getheader(imlist[0])['FILTER1']
        if filtname[0] == 'C': filtname = pyfits.getheader(imlist[0])['FILTER2']
    else: raise Exception('Instrument '+instrum+' not covered in our case list.')

    # -- init AD parameters specific to each instrument
    if instrum == 'WFC3':
        sepbits = 2048+256+64
        cr_snr = '5.0 4.0'
        skystat = 'mode'
    elif instrum == 'ACS':
        sepbits = 2048+256+64+32
        cr_snr = '4.0 3.5'
        skystat = 'mean'
        
    # -- init AD parameters based on no. of flt images
    numflts = len(imlist)
    if numflts < 4:
        cnhigh = 0
        cnlow = 0
        combtype = 'minmed'
    elif numflts < 6:
        cnhigh = 1
        cnlow = 0
        combtype = 'median'
    else:
        cnhigh = 2
        cnlow = 1
        combtype = 'median'


    # -- run AstroDrizzle
    teal.unlearn('astrodrizzle')
    if outname == 'NA': outname = filtname.lower()+'_'+finalwht.lower()
    astrodrizzle.AstroDrizzle('@imlist.dat',output=outname,num_cores=NCORES, \
                                clean=False,preserve=False, \
                                combine_type=combtype,combine_nhigh=cnhigh,combine_nlow=cnlow,
                                skywidth=0.1,skystat=skystat, \
                                driz_cr_corr=True,driz_sep_bits=sepbits,driz_cr_snr=cr_snr, \
                                final_wcs=True,final_wht_type=finalwht,final_bits=finalbits,final_rot=0.0, \
                                final_ra=finalra,final_dec=finaldec,final_outnx=finalnx,final_outny=finalny, \
                                final_scale=pixscale,final_pixfrac=pixfrac,final_kernel=finalkernel)


    # -- remove unwanted astrodrizzle files
    tmp=np.concatenate((glob.glob('*single_sci*fits'),glob.glob('*single_wht*fits'),glob.glob('*mask.fits'), \
			glob.glob('tmp*.fits'), glob.glob('*blt.fits'),glob.glob('*med.fits')))
    for ff in tmp: os.remove(ff)
