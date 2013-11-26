#! /usr/bin/env python

'''
ABOUT:
This script runs SExtractor on list of images (crclean in current directory is default).\
A trimmed catalog is constructed consisting of "max" brightest objects w/o bad flags.
The trimmed catalog is used to align flt images.

DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER for STScI, 2013

HISTORY:
Jul 2013: Original script (v0.1).
Aug 2013: Modified to run on any instrument's image (e.g., double-chip ACS/optical or single-chip WFC/IR).


FUTURE IMPROVEMENTS:


USE:
python run_sex_crclean.py                       --> make SExtractor catalogs for all crclean mages in cwd.
python run_sex_crclean.py -th 5                 --> same as above but using a 5-sigma threshold (10 is default).
python run_sex_crclean.py -max 500              --> same as above but trimming from 500 brightest objects (100 is default).
python run_sex_crclean.py -im temp.crclean.fits --> make SExtractor catalog for crclean image "tmp.crclean.fits."

'''

__author__='D.M. HAMMER'
__version__= 0.1


import os, glob, argparse, pdb, pylab, pyfits
import numpy as np
from subprocess import call


def get_wfc3_zeropoint(filter):
    # array of WFC3 zeropoints (not all here - add as needed)
    zp = {'F225W':24.0403, 'F275W':24.1305, 'F336W':24.6682, 'F390W':25.3562, 'F438W':24.8206, 'F475W':25.6799, \
	  'F547M':24.7467, 'F555W':25.7906, 'F606W':26.0691, 'F656N':20.4827, 'F658N':21.0287, 'F814W':25.0985, \
	  'F850LP':23.8338,'F105W':26.2687, 'F125W':26.2303, 'F140W':26.4524, 'F160W':25.9463}
    if zp.has_key(filter.upper()): return zp[filter.upper()]
    else: raise Exception('Zeropoint is not specified for this filter: '+filter)

def get_acs_zeropoint(hdr):
    PHOTFLAM = hdr['PHOTFLAM']
    PHOTPLAM = hdr['PHOTPLAM']
    zeropt = -2.5*np.log10(PHOTFLAM) - 5.0*np.log10(PHOTPLAM) - 2.408
    return zeropt


def make_trimmed_catalog(fullcat,trimcat,maxobj,naxis1,naxis2):
    #  trim full SEx catalog of undesireable sources for flt alignment:
    #       (a) remove objects with bad SEx flags (29 keeps 0 & 2).
    #       (b) remove objects within 32 pixels of image edge.
    #       (c) keep "maxobj" brightest objects that satisfy above criteria.

    # -- read in full catalog
    id,x,y,a,b,theta,kron,xwin,ywin,ra,dec,flux,fluxerr,mag,magerr,faper1,faper2,faper1err,faper2err,fwhm, \
       flag,sclass,isoflags,nisoflags,isoarea = np.loadtxt(fullcat,unpack=True)

    # -- select objects that satisfy criteria
    keep = np.where(((np.int32(flag) & 29) == 0) & (x>32) & (x<(naxis1-32)) & (y>32) & (y<(naxis2-32)))[0]
    fluxmin = sorted(faper2[keep],reverse=True)[np.min([len(faper2[keep]),maxobj])-1]
    gd = np.where(((np.int32(flag) & 29) == 0) & (x>32) & (x<(naxis1-32)) & (y>32) & (y<(naxis2-32)) & \
                 (faper2 >= fluxmin))[0]

    # -- save trimmed catalog
    np.savetxt(trimcat,zip(id[gd],x[gd],y[gd],a[gd],b[gd],theta[gd],kron[gd],xwin[gd],ywin[gd],ra[gd],dec[gd], \
       flux[gd],fluxerr[gd],mag[gd],magerr[gd],faper1[gd],faper2[gd],faper1err[gd],faper2err[gd],flag[gd], \
       sclass[gd],isoflags[gd],nisoflags[gd],isoarea[gd]))

    return len(gd)



if __name__=='__main__':

    # -- Parse input parameters
    #--------------------------
    parser = argparse.ArgumentParser(description='Run SExtractor on specified images.')
    parser.add_argument('-im', '--images',default='*crclean.fits', type=str, help='Input fits file(s). \
                         Default is all CR-cleaned science images in working directory.')
    parser.add_argument('-max', '--maxobj',default=100, type=int, help='Input maximum number of SExtractor detections \
                         to use for tweakreg matching (each chip). Default is 100 objects per chip.')
    parser.add_argument('-th', '--threshold',default=10, type=float, help='Input threshold used for BOTH detect_thresh \
                         and analysis_thresh. Default is 10 as tests showed better tweakreg residuals for high thresh.')

    options = parser.parse_args()
    imlist = glob.glob(options.images)
    imlist.sort()
    maxobj = options.maxobj
    threshold = options.threshold


    # -- initialize variables to hold names of catalogs for each image
    tcatnamesA_global = []
    tcatnamesB_global = []


    # -- Run SExtractor on each image
    #--------------------------------
    for im in imlist:

        # -- get filter name, exposure time, zeropoint, & gain value
        fheader = pyfits.getheader(im)
        instr = fheader['INSTRUME']
        detector = fheader['DETECTOR']
        if instr == 'WFC3':
	    #  ***NOTE: WFC3 IR crclean images are in cnts, while flts are in cnts/sec
            filtname = fheader['FILTER']
	    exptime = fheader['EXPTIME']
            zeropt = get_wfc3_zeropoint(filtname) + 2.5*np.log10(exptime)
            hstgain = 1.0
            if detector == 'UVIS': pscale = 0.03962
            elif detector == 'IR': pscale = 0.12825
            else: raise Exception('Detector '+detector+' not covered in our case list.')
            naxis1 = pyfits.getval(im,'NAXIS1',ext=('SCI',1))
            naxis2 = pyfits.getval(im,'NAXIS2',ext=('SCI',1))
            configname = instr.lower()+detector.upper()+'.sex.crclean.config'
        elif instr == 'ACS':
            filtname = fheader['FILTER1']
            if filtname[0] == 'C': filtname = fheader['FILTER2']
	    exptime = fheader['EXPTIME']
	    scihdr = pyfits.getheader(im,ext=('SCI',1))
            zeropt = get_acs_zeropoint(scihdr) + 2.5*np.log10(exptime)
            hstgain = fheader['CCDGAIN']
            if detector == 'WFC': pscale = 0.049
            else: raise Exception('Detector '+detector+' not covered in our case list.')
            naxis1 = pyfits.getval(im,'NAXIS1',ext=('SCI',1))
            naxis2 = pyfits.getval(im,'NAXIS2',ext=('SCI',1))
            configname = instr.lower()+detector.upper()+'.sex.crclean.config'
        else: raise Exception('Instrument '+instr+' is not supported in case list.')


	# -- determine the science extension, store chip information, & assign output catalog name
	sciext = []
	chipid = []
	catname = []
	tcatname = []
        hdulist = pyfits.open(im)
        for ff in xrange(len(hdulist)):
            if hdulist[ff].name == 'SCI':
                sciext.append(ff)
                fheader = pyfits.getheader(im,ext=sciext[-1])
		chipid.append(fheader.get('CCDCHIP',default=1))
		catname.append(im.split('.fits')[0] +'_sci'+str(len(chipid))+'.sex.all')
                tcatname.append(catname[-1].split('.all')[0])
		if len(sciext) == 1: tcatnamesA_global.append(tcatname[-1])
		elif len(sciext) == 2: tcatnamesB_global.append(tcatname[-1])
		else: raise Exception('Unexpected number of SCI extensions (>2).')
	if len(sciext) == 1: tcatnamesB_global.append('')


        # -- run SExtractor
	for ff in xrange(len(sciext)):

	    # -- create SE catalogs & make corresponding ds9 region file with Kron apertures.
	    #     NOTE: SExtractor considers the zeroth extension to be first fits ext with image data.
	    call(['sex','-c',configname,im+'['+str(sciext[ff]-1)+']','-CATALOG_NAME',catname[ff], \
	          '-PIXEL_SCALE',str(pscale),'-GAIN',str(hstgain),'-MAG_ZEROPOINT',str(zeropt), \
	          '-DETECT_THRESH',str(threshold),'-ANALYSIS_THRESH',str(threshold),'-FLAG_IMAGE','','-FLAG_TYPE',''])
	    call(['cat2reg',catname[ff]])


	    # -- trim catalog of undesireable sources for flt alignment; make ds9 file for trimmed catalog
            nobj = make_trimmed_catalog(catname[ff],tcatname[ff],maxobj,naxis1,naxis2)
            call(['cat2reg',tcatname[ff]])
            print '# OBJECTS IN TRIMMED CATALOG: '+str(nobj)



    # -- Create "catfile" for input to tweakreg
    #----------------------------------------
    fltlist = [imlist[ff].split('crclean')[0] + 'flt.fits' for ff in xrange(len(imlist))]
    np.savetxt('catfile.sex',zip(fltlist,tcatnamesA_global,tcatnamesB_global),fmt='%s')

