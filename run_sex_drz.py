#! /usr/bin/env python

'''
ABOUT:
This script runs SExtractor on a list of drizzled images (no default - must input filenames).
A trimmed catalog is constructed consisting of "maxobj" brightest objects w/o bad flags, with
the option of also selecting on fwhm_image and class_star.
The trimmed catalog is used to co-align drizzled images.

DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER, 2013

HISTORY:
Sep 2013: Original script (v0.1).
Oct 2013: Changed default RMS value for zero-exposure pixels to 1e32.

FUTURE IMPROVEMENTS:


USE:
python run_sex_drz.py -im drz.fits                         -->construct SEx catalogs on drizzled image "drz.fits."
python run_sex_drz.py -im drz.fits -ivm wht.fits           -->use IVM wht map (converted to rms map) for weighting.
python run_sex_drz.py -im drz.fits -exp exp.fits           -->use exposure map for weighting and for exposure flag map.
python run_sex_drz.py -im drz.fits -exp e.fits -ivm i.fits -->use IVM for SEx weight & make exposure flag map.
python run_sex_drz.py -im drz.fits -max 500 -th 3          -->trim catalog using 500 brightest obj; use Sex thresh=3.
python run_sex_drz.py -im drz.fits -cs 0.90 -f 2.0         -->trim catalog using w/CLASS_STAR>0.9 && FWHM_IMAGE<2.


'''

__author__='D.M. HAMMER'
__version__= 0.2


import os, glob, argparse, pdb, pylab, pyfits, shutil
import numpy as np
from subprocess import call
from scipy.ndimage.filters import median_filter
from scipy.ndimage.filters import uniform_filter


def get_wfc3_zeropoint(filter):
    # array of WFC3 AB zeropoints (not all here - add as needed)
    zp = {'F225W':24.0403, 'F275W':24.1305, 'F336W':24.6682, 'F390W':25.3562, 'F438W':24.8206, 'F475W':25.6799, \
          'F547M':24.7467, 'F555W':25.7906, 'F606W':26.0691, 'F656N':20.4827, 'F658N':21.0287, 'F814W':25.0985, \
          'F850LP':23.8338,'F105W':26.2687, 'F125W':26.2303, 'F140W':26.4524, 'F160W':25.9463}

    if zp.has_key(filter.upper()): return zp[filter.upper()]
    else: raise Exception('Zeropoint is not specified for this filter: '+filter)

def get_acs_zeropoint(hdr):
    PHOTFLAM = hdr['PHOTFLAM']
    PHOTPLAM = hdr['PHOTPLAM']
    zeropt = -2.5*np.log10(PHOTFLAM) - 5.0*np.log10(PHOTPLAM) - 2.408  #AB
    return zeropt


def make_rmsmap(whtim, filtname):
    rmsname = filtname.lower()+'_rms.fits'
    shutil.copy(whtim,'./'+rmsname)
    rms = pyfits.open(rmsname,mode='update')
    data = rms[0].data
    rms[0].data = 1.0/np.sqrt(data)			   # assumes whtim is inverse-variance map
    rms[0].data[np.isfinite(rms[0].data) == False] = 1e32  # assign large number above default WEIGHT_THRESH (1e30)
    rms.close()
    return rmsname


def make_flagmap(expim,filtername):
    ''' make flag map based on [WHT] exposure map: 0=good,1=low exposure (<2/3 median),2=edge,4=zero exposure'''
    zeroexpflag = 4
    edgegrowflag = 2
    lowexpflag = 1
    exp = pyfits.getdata(expim)

    # -- construct flag map just for border pixels
    edgemap = make_edge_flagmap(exp,zflag=zeroexpflag,gflag=edgegrowflag)

    # -- construct flag map just for interior low-exposure pixels
    lowexpmap = make_lowexp_flagmap(exp,0.6667*np.median(exp[exp > 0]),lflag=lowexpflag)

    # -- construct final consolidated flag map
    flagmap = edgemap + lowexpmap
    flagmap[edgemap == zeroexpflag] = zeroexpflag	#re-set zero-exposure pixels

    # -- save flag image
    flagname = filtername.lower()+'_flagmap.fits'
    shutil.copy(expim,'./'+flagname)
    flaghdu = pyfits.open(flagname,mode='update')
    flaghdu[0].data = flagmap.astype('int32')
    flaghdu.close()
    return flagname			#return name of flagmap

def make_edge_flagmap(exp,thresh=0.0,boxsz=67,zflag=4,gflag=2):
    # -- identify pixels with zero exposure (less than thresh)
    zeromap = np.zeros(exp.shape)
    zeromap[exp <= thresh] = 1.0

    # -- median filter zero-exposure map to reject pixels that have little exposure owing to heavily-cleaned CRs
    medf_zeromap = median_filter(zeromap,size=15,mode='constant',cval=1.0)

    # -- boxcar smooth the median-filtered image to "grow" image edge
    smooth_medfmap = uniform_filter(medf_zeromap,size=boxsz,mode='constant',cval=1.0)
    smooth_medfmap[smooth_medfmap < 1e-06] = 0.0

    # -- combine zero and near-zero edge pixels into single flag map
    edge_flagmap = np.zeros(exp.shape)
    edge_flagmap[exp <= thresh] = zflag
    edge_flagmap[exp > thresh] = smooth_medfmap[exp > thresh]
    edge_flagmap[(edge_flagmap > 0.0) & (edge_flagmap < 1.0)] = gflag
    return edge_flagmap

def make_lowexp_flagmap(exp,thresh,boxsz=7,lflag=1):
    # -- median filter image to reject localized low-exposure pixels owing to heavily-cleaned CRs
    medexp = median_filter(exp,size=boxsz,mode='constant',cval=1.0)

    # -- convert median-filtered map to simple flag image (0>thresh; 1<thresh)
    smedexp = np.zeros(exp.shape)
    smedexp[medexp <= thresh] = 1.0

    # -- smooth simple flag image and identify non-zero pixels (i.e., low-exposure pixels)
    boxsm_smedexp = uniform_filter(smedexp,size=boxsz,mode='constant',cval=1.0)
    boxsm_smedexp[boxsm_smedexp < 1e-06] = 0.0

    # -- build low-exposure flagmap
    lowexp = np.zeros(exp.shape)
    lowexp[boxsm_smedexp > 0.0] = lflag
    return lowexp


def make_trimmed_catalog(fullcat,trimcat,maxobj,maxfwhm,minstellar):
    #  trim full SEx catalog of undesireable sources for drz to drz alignment:
    #       (a) remove objects with bad SEx flags (29 keeps 0 & 2).
    #       (b) remove objects with bad exposure flags (keep 0; 1 only if <10% pixels are affected).
    #       (c) optionally select on size/profile via fwhm and/or class_star.
    #       (d) keep "maxobj" brightest objects that satisfy above criteria.

    # -- read in full catalog
    id,x,y,a,b,theta,kron,xwin,ywin,ra,dec,flux,fluxerr,mag,magerr,faper1,faper2,faper1err,faper2err,fwhm, \
       flag,sclass,isoflags,nisoflags,isoarea = np.loadtxt(fullcat,unpack=True)

    # -- select objects that satisfy criteria
    keep = np.where(((np.int32(flag) & 29) == 0) & (isoflags < 2) &   \
                    (nisoflags.astype('float32')/isoarea.astype('float32') <= 0.1) & \
                    (sclass >= minstellar) & (fwhm <= maxfwhm))[0]
    fluxmin = sorted(faper2[keep],reverse=True)[np.min([len(faper2[keep]),maxobj])-1]
    gd = np.where(((np.int32(flag) & 29) == 0) & (isoflags < 2) & \
                  (nisoflags.astype('float32')/isoarea.astype('float32') <= 0.1) & \
                  (sclass >= minstellar) & (fwhm <= maxfwhm) & (faper2 >= fluxmin))[0]

    # -- save trimmed catalog
    np.savetxt(trimcat,zip(id[gd],x[gd],y[gd],a[gd],b[gd],theta[gd],kron[gd],xwin[gd],ywin[gd],ra[gd],dec[gd], \
       flux[gd],fluxerr[gd],mag[gd],magerr[gd],faper1[gd],faper2[gd],faper1err[gd],faper2err[gd],flag[gd], \
       sclass[gd],isoflags[gd],nisoflags[gd],isoarea[gd]))
    return len(gd)


if __name__=='__main__':

    # -- Parse input parameters
    #--------------------------
    parser = argparse.ArgumentParser(description='Run SExtractor on specified images.')
    parser.add_argument('-im', '--images',default='NONE', type=str, help='Input drizzled fits image(s) for photometry. \
                         Default is NONE - must input image(s).')
    parser.add_argument('-dim', '--detectim',default='NONE', type=str, help='Input drizzled fits image(s) for detection. \
                         Default is NONE - run SExtractor using same input image for both detection & photometry.') 
    parser.add_argument('-ivm', '--ivmname',default='NONE', type=str, help='Input IVM weight map image(s). \
                         Default is NONE. Number of weight maps must match number of drizzled images.')
    parser.add_argument('-exp', '--expname',default='NONE', type=str, help='Input effective exposure map. \
                         Default is NONE. Number of exposure maps must match number of drizzled images. \
                         IVM switch takes precedence, i.e., if both -ivm & -exp, will drizzle with IVM final_wht_type \
                         but will use exposure map to construct flag map to exclude bad sources in trimmed catalog.')
    parser.add_argument('-th', '--threshold',default=10, type=float, help='Input threshold used for BOTH detect_thresh \
                         and analysis_thresh. Default is 10 - tests showed better tweakreg residuals at 10.')
    parser.add_argument('-max', '--maxobj',default=300, type=int, help='Input maximum number of SExtractor detections \
                         to use for tweakreg matching. Default is 300 objects per extension.')
    parser.add_argument('-cs', '--stellarity',default=0.0, type=float, help='Input lower limit for selecting objects \
                         via SExtractor "CLASS_STAR." Default is 0.0, i.e., do not select on stellarity.')
    parser.add_argument('-f', '--fwhm',default=100.0, type=float, help='Input upper limit for selecting objects via \
                         SExtractor "FWHM_IMAGE." Default is 100, i.e., do not select on FWHM.')
    options = parser.parse_args()
    threshold = options.threshold
    maxobj = options.maxobj
    minstellar = options.stellarity
    maxfwhm = options.fwhm   

    # -- initialize input image list
    if options.images == 'NONE': raise Exception('Must input drizzled images.')
    imlist = glob.glob(options.images)
    imlist.sort()
    
    if options.detectim != 'NONE':
        dimlist = glob.glob(options.detectim)
        dimlist.sort()
        # case of many input images and only one detection image
        if ((len(dimlist) == 1) & (len(dimlist) < len(imlist))):
            tmpdim = dimlist[0]
            dimlist = []
            for ff in imlist: dimlist.append(tmpdim)
    else:
        # use same image for detection and photometry
        dimlist = []
        for ff in imlist: dimlist.append(ff)	


   # -- initialize weight,exposure image lists & weighting type
    if options.ivmname != 'NONE':
        whtlist = glob.glob(options.ivmname)
        whtlist.sort()
	if len(whtlist) != len(imlist): raise Exception('Number of ivm maps must equal number of drizzled images.')
	weighttype ='MAP_RMS'

    if options.expname != 'NONE':
        explist = glob.glob(options.expname)
        explist.sort()
        if len(explist) != len(imlist): raise Exception('Number of exposure maps must equal number of drizzled images.')
        if options.ivmname == 'NONE':
            whtlist = explist
	    weighttype ='MAP_WEIGHT'
    else:
        explist = []
        for ff in imlist: explist.append('NONE')

    if ((options.ivmname == 'NONE') & (options.expname == 'NONE')):
        weighttype ='NONE'
	whtlist = []
	for ff in imlist: whtlist.append('NONE')


    tcatname_all = []
  
    # -- Run SExtractor on each image
    #--------------------------------
    for im,dim,wht,exp in zip(imlist,dimlist,whtlist,explist):

	# -- assign various image properties (gain,exposure,zeropt,filter name, pixel scale)
        fheader = pyfits.getheader(im)
        instr = fheader['INSTRUME']
        if instr == 'WFC3':
            filtname = fheader['FILTER']
	    exptime = fheader['EXPTIME']
	    zeropt = get_wfc3_zeropoint(filtname)
	    hstgain = 1.0
	    gain = hstgain * exptime
	    pscale = fheader['D001SCAL']
	    detector = pyfits.getval(im,'DETECTOR',ext=0)
            configname = instr.lower()+detector.upper()+'.sex.drz.config'
            if detector == 'IR': saturation = 10000000.    # assign large value - IR images aren't saturated
            elif detector == 'UVIS':
                # use median exposure time to estimate saturation level
	        ftable = pyfits.getdata(im,ext=1)
	        medexp = np.median(ftable['EXPTIME'])
	        saturation = 80000./medexp
	    else: raise Exception('Detector '+detector+' is not covered in our case list.')
        elif instr == 'ACS':
            filtname = fheader['FILTER1']
            if filtname[0] == 'C': filtname = fheader['FILTER2']
            exptime = fheader['EXPTIME']
	    zeropt = get_acs_zeropoint(fheader)
	    hstgain = fheader['CCDGAIN']
	    gain = hstgain * exptime
	    # use median exposure time to estimate saturation level
	    ftable = pyfits.getdata(im,ext=1)
	    medexp = np.median(ftable['EXPTIME'])
	    saturation = 80000./medexp
	    pscale = fheader['D001SCAL']
	    detector = pyfits.getval(im,'DETECTOR',ext=0)
            configname = instr.lower()+detector.upper()+'.sex.drz.config'
	else: raise Exception('Instrument '+instr+' is not covered in our case list.')


	# -- assign output catalog name
	catname = im.split('fits')[0]+'sex.all'
	tcatname = im.split('fits')[0]+'sex'
	tcatname_all.append(tcatname)


	# -- construct rms map if "IVM" weight map was input.
	if weighttype == 'MAP_RMS':
            rmsname = make_rmsmap(wht,filtname)
            wht = rmsname


	# -- construct flag map if exposure map was input.
	if options.expname != 'NONE':
	    flagname = make_flagmap(exp,filtname)
	    flagtype = 'OR'
	else:
	    flagname = ''
	    flagtype = ''	    


        # -- call SExtractor & make ds9 region file for full catalog
        call(['sex',dim+','+im,'-c',configname,'-WEIGHT_IMAGE',wht+','+wht,'-WEIGHT_TYPE',weighttype+','+weighttype, \
              '-CATALOG_NAME',catname,'-PIXEL_SCALE',str(pscale),'-GAIN',str(gain),'-MAG_ZEROPOINT',str(zeropt), \
              '-SATUR_LEVEL',str(saturation),'-DETECT_THRESH',str(threshold),'-ANALYSIS_THRESH',str(threshold), \
              '-FLAG_IMAGE',flagname,'-FLAG_TYPE',flagtype])
	call(['cat2reg',catname])


	# -- trim catalog of undesireable sources for drz to drz alignment; make ds9 region file for trimmed catalog
        nobj = make_trimmed_catalog(catname,tcatname,maxobj,maxfwhm,minstellar)
        print '# OBJECTS IN TRIMMED CATALOG: '+str(nobj)
        call(['cat2reg',tcatname])
       	

    # -- Create "catfile" for input to run_tweakreg_drz
    np.savetxt('catfile.sex',zip(imlist,tcatname_all),fmt='%s')

