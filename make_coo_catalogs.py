#! /usr/bin/env python

'''
ABOUT:
This program makes star catalogs from crclean.fits images using the tweakreg task in Drizzlepac.

DEPENDS:
os
glob
pyfits
drizzlepac

AUTHOR:
Roberto J. Avilas

HISTORY:
August 2013: Original script (v0.1).

FUTURE IMPROVEMENTS:
Need to output catalogs in order of brightness

USE:
python make_catalogs.py

'''

__author__='Roberto J. Avila'
__version__= 0.1


import glob
import os
import pyfits
import numpy as np
from drizzlepac import tweakreg as tweak

imlist = glob.glob('*crclean.fits')
suffix1 = 'sci1.coo'
suffix2 = 'sci2.coo'

for i,im in enumerate(imlist):

    hdu = pyfits.open(im)
    hdr0 = hdu[0].header
    hdr1 = hdu[1].header
    hdu.close

    camera = hdr0['detector']

    if camera == 'WFC':
        conv_w = 3.5
        pkmax = 70000.
        pkmin = 500. 
        
    elif camera == 'IR':
        conv_w = 2.5
        pkmax = 70000.
        pkmin = 1000.   
       
    else:
        print 'Cannot figure out camera for {}'.format(im)
        continue


    root = im[:10]
    tweak.TweakReg(im,\
                   writecat=True,\
                   clean=True,\
                   runfile='test.log',\
                   conv_width=conv_w,\
                   updatehdr=False,\
                   headerlet=False,\
                   shiftfile=False,\
                   catfile='',\
                   computesig=True,\
                   peakmax=pkmax,\
                   peakmin=pkmin,\
                   threshold=10.)


    columns = hdr1['naxis1']
    rows = hdr1['naxis2']
    redge = columns - 50
    ledge = 50
    tedge = rows - 50
    bedge = 50

    x1,y1 = np.loadtxt(root+'crclean_sci1_xy_catalog.coo',\
        usecols=(0,1),unpack = True)
    ind1 = [(x1>ledge)&(x1<redge)&(y1>bedge)&(y1<tedge)]
    np.savetxt(root+suffix1,zip(x1[ind1],y1[ind1]))

    os.remove(root+'crclean_refxy_catalog.coo')
    os.remove(root+'crclean_sky_catalog.coo')
    os.remove(root+'crclean_sci1_xy_catalog.coo')

    if camera == 'WFC':
        x2,y2 = np.loadtxt(root+'crclean_sci2_xy_catalog.coo',\
            usecols=(0,1),unpack = True)
        ind2 = [(x2>ledge)&(x2<redge)&(y2>bedge)&(y2<tedge)]
        np.savetxt(root+suffix2,zip(x2[ind2],y2[ind2]))
        os.remove(root+'crclean_sci2_xy_catalog.coo')


print 'Writing catfile now\n'
catfile = open('catfile.coo','w')

for i in range(len(imlist)):
    root = imlist[i][:10]
    if camera == 'WFC':catfile.write('{0} {1} {2}\n'.format(root+'flc.fits',root+suffix1,root+suffix2))
    if camera == 'IR' :catfile.write('{0} {1}\n'.format(root+'flt.fits',root+suffix1))

catfile.close()

