#! /usr/bin/env python

'''
ABOUT:
This script reads in the matched catalogs from tweakreg and recreates the four
diagnostic delta-x/y vs. x/y residual diagrams, and the single vectorgram, constructed by tweakreg in interactive mode.
This is necessary as tweakreg does not allow us to save the diagnostic diagrams in batch mode.


DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER, 2013

HISTORY:
July 2013: Original script (v0.1).
Oct. 2013: Added Roberto fix so matplotlib is run in batch mode (previously interactive which failed during multithread).
Oct. 2013: Removed astropy import (not required). Added new "rootname" variable to label the source image on diagram.
           Changed the method of selecting the output name for diagrams (image name minus .fits + "_xyresids.pdf").


FUTURE IMPROVEMENTS:
Correct title information for files with LONG names.
Accept reference image to get proper scaling information (or at least a name of ref catalog).


USE:
python run_tweakreg_make_xyresids.py                            -->constructs diagrams for each *fit.match file in cwd.
python run_tweakreg_make_xyresids.py -cat test1.fit.match       -->same as above but only for file "test1.fit.match."
'''

__author__='D.M. HAMMER'
__version__= 0.1


import glob, argparse, pdb, pyfits
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.ticker import MultipleLocator
import pylab


def parse_tweakcat(catname):
    f = open(catname)
    xref = []
    yref = []
    xs = []
    ys = []

    for line in f:
        if line[0] != '#':
	    xref.append(line[0:15])
	    yref.append(line[15:30])
	    xs.append(line[90:105])
	    ys.append(line[105:120])

    xref = np.array(xref).astype('float32')
    yref = np.array(yref).astype('float32')
    xs = np.array(xs).astype('float32')
    ys = np.array(ys).astype('float32')
    f.close()
    return xref,yref,xs,ys


if __name__=='__main__':
    # -- parse input parameters
    parser = argparse.ArgumentParser(description='Read in xy residuals from tweakreg matched catalog files and plot.')
    parser.add_argument('-cat','--catalogs',default='*fit.match', type=str, help='Input tweakreg catalogs with residuals.\
                        Default is all *fit.match files in working directory.')
    options = parser.parse_args()


    # -- save images and tweakreg catalogs in arrays
    fitcats = glob.glob(options.catalogs)
    fitcats.sort()
    images = [fitcats[i].split('_catalog')[0]+'.fits' for i in xrange(len(fitcats))]
    impref = [images[i].split('.fits')[0] for i in xrange(len(images))]
    rootname = [pyfits.getval(images[i],'ROOTNAME',ext=0).strip() for i in xrange(len(images))]


    # Iterate over each catalog and make diagrams
    # --------------------------------------------
    for cat,im,pref,root,ctr in zip(fitcats,images,impref,rootname,xrange(len(fitcats))):
    	
        # -- record image properties if available (covers both drz and flt)
        checkim = glob.glob(im)
        if len(checkim) > 0:
            fheader = pyfits.getheader(im)
            instr = fheader['INSTRUME']
            detector = fheader['DETECTOR']

            # -- record filter name
            if instr == 'WFC3':
                filtname = fheader['FILTER']
            elif instr == 'ACS':
                filtname = fheader['FILTER1']
                if filtname[0] == 'C': filtname = fheader['FILTER2']
            else: raise Exception('Instrument '+instr+' not covered in our case list.')

            # -- record pixel scale and image size
            pscale = fheader.get('D001SCAL',default='NA')
            if pscale != 'NA':
                imtype = 'drz'
	        naxis1 = fheader['NAXIS1']
	        naxis2 = fheader['NAXIS2']
            else:
		imtype ='flt'
		if detector == 'UVIS': pscale = 0.03962
		elif detector == 'IR': pscale = 0.12825
		elif detector == 'WFC': pscale = 0.049
		else: raise Exception('Detector '+detector+' not covered in our case list.')
		naxis1 = pyfits.getval(im, 'NAXIS1', ext=('SCI',1))
		naxis2 = pyfits.getval(im, 'NAXIS2', ext=('SCI',1))
	else: raise Exception('An image with the root name '+im+' must be located in cwd.')


        # -- read in matched catalog (note that residuals are fit-input)
	xref,yref,xs,ys = parse_tweakcat(cat)


        # -- generate fit statistics
        xrms = np.std(xs)
        yrms = np.std(ys)
	xrms_mas = xrms * pscale * 1000.0
	yrms_mas = yrms * pscale * 1000.0


	# -- establish x-range & xtick intervals for diagrams
	if ((imtype == 'flt') & (detector == 'IR')):
	    # case of WFC3 IR flt
            x0 = -520
	    x1 = 520
	    xtmin = 50
	    xtmax = 200
            xtmin_vec = 50
            xtmax_vec = 200
            ytmin = 0.1
            ytmax = 0.2
	elif imtype == 'flt':
	    # case of ACS/UVIS flt
            x0 = -2100
            x1 = 2100
	    xtmin = 100
	    xtmax = 1000
            xtmin_vec = 200
            xtmax_vec = 1000
            ytmin = 0.1
            ytmax = 0.2
	elif imtype == 'drz':
            axsize_min = (np.max([naxis1,naxis2]) * 1.025 + 20.0 - np.remainder(np.max([naxis1,naxis2]) * 1.025,20))/2.0
            axsize_max = (np.max([np.max(np.abs(xref)),np.max(np.abs(yref))]) * 1.025 + 20.0 - \
                          np.remainder(np.max([np.max(np.abs(xref)),np.max(np.abs(yref))]) * 1.025,20))
            axsize = np.max([axsize_min,axsize_max])
            x0 = -1.0*axsize
            x1 = axsize
            xtmin = 100
	    xtmax = 1000
            xtmin_vec = 100
            xtmax_vec = 500
            ytmin = 0.2
            ytmax = 0.4

        # -- establish y-range for diagrams
        ymax = np.max(np.abs([np.min(np.concatenate((xs,ys))),np.max(np.concatenate((xs,ys)))]))
        y0 = -1.0 * np.max([(ymax + 0.1),0.5])
        y1 = np.max([(ymax + 0.1),0.5])


	# MAKE X/Y RESIDUALS VS X/Y POSITION DIAGRAM
        # ------------------------------------------
	if ctr == 0:
            fig1=pylab.figure(figsize=(11,8))
            fig1.subplots_adjust(wspace=0.3, hspace=0.2)
            fig1.suptitle(filtname+' XY Residuals ('+root+').  RMS(X)={:1.2f} pix ({:2.2f} mas)'.format(xrms,xrms_mas)+ \
                          '  RMS(Y)={:1.2f} pix ({:2.2f} mas)'.format(yrms,yrms_mas)+'  # objects='+str(len(xs)), \
                          ha='center', color='black', weight='normal')
	else:
	    pylab.figure(fig1.number)
	    pylab.clf()
	    fig1.suptitle(filtname+' XY Residuals ('+root+').  RMS(X)={:1.2f} pix ({:2.2f} mas)'.format(xrms,xrms_mas)+ \
                          '  RMS(Y)={:1.2f} pix ({:2.2f} mas)'.format(yrms,yrms_mas)+'  # objects='+str(len(xs)), \
                          ha='center', color='black', weight='normal')

	# -- delta-X vs. X
        ax1=pylab.subplot(2,2,1)
        ax1.yaxis.set_minor_locator(MultipleLocator(ytmin))
        ax1.yaxis.set_major_locator(MultipleLocator(ytmax))
        ax1.xaxis.set_minor_locator(MultipleLocator(xtmin))
        ax1.xaxis.set_major_locator(MultipleLocator(xtmax))

	pylab.scatter(xref,xs,s=4)
	pylab.axhline(0.0,color='red')
	pylab.xlim(x0,x1)
	pylab.ylim(y0,y1)
	pylab.xlabel('X (pixels)')
	pylab.ylabel('$\Delta$X (pixels)')


        # -- delta-X vs. Y
        ax1=pylab.subplot(2,2,2)
        ax1.yaxis.set_minor_locator(MultipleLocator(ytmin))
        ax1.yaxis.set_major_locator(MultipleLocator(ytmax))
        ax1.xaxis.set_minor_locator(MultipleLocator(xtmin))
        ax1.xaxis.set_major_locator(MultipleLocator(xtmax))

        pylab.scatter(yref,xs,s=4)
        pylab.axhline(0.0,color='red')
        pylab.xlim(x0,x1)
        pylab.ylim(y0,y1)
        pylab.xlabel('Y (pixels)')
        pylab.ylabel('$\Delta$X (pixels)')


        # -- delta-Y vs. X
        ax1=pylab.subplot(2,2,3)
        ax1.yaxis.set_minor_locator(MultipleLocator(ytmin))
        ax1.yaxis.set_major_locator(MultipleLocator(ytmax))
        ax1.xaxis.set_minor_locator(MultipleLocator(xtmin))
        ax1.xaxis.set_major_locator(MultipleLocator(xtmax))
        
        pylab.scatter(xref,ys,s=4)
        pylab.axhline(0.0,color='red')
        pylab.xlim(x0,x1)
        pylab.ylim(y0,y1)
        pylab.xlabel('X (pixels)')
        pylab.ylabel('$\Delta$Y (pixels)')


        # -- delta-Y vs. Y
        ax1=pylab.subplot(2,2,4)
        ax1.yaxis.set_minor_locator(MultipleLocator(ytmin))
        ax1.yaxis.set_major_locator(MultipleLocator(ytmax))
        ax1.xaxis.set_minor_locator(MultipleLocator(xtmin))
        ax1.xaxis.set_major_locator(MultipleLocator(xtmax))

        pylab.scatter(yref,ys,s=4)
        pylab.axhline(0.0,color='red')
        pylab.xlim(x0,x1)
        pylab.ylim(y0,y1)
        pylab.xlabel('Y (pixels)')
        pylab.ylabel('$\Delta$Y (pixels)')

	# -- save figure
	pylab.savefig(pref+'_xyresids.pdf')
	pylab.savefig(pref+'_xyresids.png')



        # MAKE VECTORGRAM
	# ----------------
        if ctr == 0: fig2=pylab.figure()
	else:
            pylab.figure(fig2.number)
            pylab.clf()

	# -- determine the best scale to use (must adjust for large offsets, e.g., ACS-->WFC3, or vectors are huge).
	#if np.max([yrms,xrms]) > 0.2:
        if imtype == 'drz':
            vscale = 1.0/500.0
            legname = '0.2 pix'
            legsize = 0.2    
        else:
            vscale = 1.0/2000.0
            legname = '0.05 pix'
            legsize = 0.05

        # -- plot vector diagram
        ax2=pylab.subplot(1,1,1,aspect='equal')
        Q=ax2.quiver(xref,yref,xs,ys,color='r', angles='xy',units='xy',scale=vscale)
        pylab.xlabel('X')
        pylab.ylabel('Y')
        pylab.xlim(x0-100,x1+100)
        pylab.ylim(x0-100,x1+100)
        ax2.yaxis.set_minor_locator(MultipleLocator(xtmin_vec))
        ax2.yaxis.set_major_locator(MultipleLocator(xtmax_vec))
        ax2.xaxis.set_minor_locator(MultipleLocator(xtmin_vec))
        ax2.xaxis.set_major_locator(MultipleLocator(xtmax_vec))
        pylab.title(filtname+' XY Residuals ('+root+').  RMS(X)={:1.2f} pix ({:2.2f} mas)'.format(xrms,xrms_mas)+ \
        	    '  RMS(Y)={:1.2f} pix ({:2.2f} mas)'.format(yrms,yrms_mas)+'  # objects='+str(len(xs)),size='small')
        qk = pylab.quiverkey(Q, 0.915, 0.937,legsize,legname, labelpos='N', coordinates='axes', \
        		     fontproperties={'weight': 'bold'}, color='black')
        pylab.show()
        pylab.savefig(pref+'_xyvector.pdf')
        pylab.savefig(pref+'_xyvector.png')

