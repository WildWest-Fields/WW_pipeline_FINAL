#! /usr/bin/env python

'''
ABOUT:
This programs constructs crclean images from a list of ASN files or flat-fielded HST images.

DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER, 2013

HISTORY:
August 2013: Original script (v0.1).
October 2013: Removed iraf unlearn (unnecessary). Changed updatewcs=True as FF uses unprocessed raws.
              Set context=False in call to AD (no need for ctx files when using crcleans).
              
FUTURE IMPROVEMENTS:

USE:
python run_make_crclean.py              --> make crcleans from flts located in asn.fits files in cwd (DEFAULT).
python run_make_crclean.py -na          --> make crcleans using all flts/flcs in current directory, in one call.

'''

__author__='D.M. HAMMER'
__version__= 0.1

import os, glob, argparse, pyfits, pdb
import numpy as np
from drizzlepac import astrodrizzle
from stsci.tools import teal


if __name__=='__main__':
    # -- Parse input parameters
    parser = argparse.ArgumentParser(description='Run tweakreg on input images using custom source catalogs.')
    parser.add_argument('-a', '--asn',default='*asn.fits', type=str, help='Input association fits file(s). \
                                Default is all _asn files in current directory.')
    parser.add_argument('-na', '--no_asn',default=False, action='store_true', help='Should we construct crcleans using \
                                all flts in current directory in one call? Default is False, \
                                i.e., by default, we construct crcleans only for flts in each _asn.fits file.')
    parser.add_argument('-im', '--images',default='*fl?.fits', type=str, help='Input fl? fits images. \
                                Default is all _fl? images in current directory.')                                
    parser.add_argument('-log', '--logfile',default='crclean.log', type=str, help='Input name of astrodrizzle log file. \
                                Default is crclean.log')
    parser.add_argument('-c', '--cores',default=4, type=int, help='Input number of processing cores used by AD. \
                                Default is 4.')

    options = parser.parse_args()
    USE_ASN = not(options.no_asn)
    NCORES = options.cores
    logfile = options.logfile


    # -- construct crclean images
    teal.unlearn('astrodrizzle')
    
    if USE_ASN:
        asnlist = glob.glob(options.asn)
        asnlist.sort()
        if len(asnlist) == 0: raise Exception('No asn files located. Use switch "-na" if using array of flt images.')
	for asn in asnlist: astrodrizzle.AstroDrizzle(asn,num_cores=NCORES,preserve=False,runfile=logfile,context=False,\
                                  updatewcs=True,driz_sep_bits='2048,256,64,32',driz_cr_corr=True,driz_combine=False)
                                  
    else:
        imlist = glob.glob(options.images)
        imlist.sort()
        f = open('imlist.dat', 'w')
        for ff in imlist: f.write(ff+'\n')
        f.close()
	astrodrizzle.AstroDrizzle('@imlist.dat',num_cores=NCORES,preserve=False,runfile=logfile,context=False,\
                                  updatewcs=True,driz_sep_bits='2048,256,64,32',driz_cr_corr=True,driz_combine=False)


    # -- remove unwanted files
    bad = np.concatenate((glob.glob('*single*fits'),glob.glob('*crmask.fits'),glob.glob('tmp*.fits'), \
                          glob.glob('*blt.fits'), glob.glob('*med.fits')))
    for tmp in bad: os.remove(tmp)
