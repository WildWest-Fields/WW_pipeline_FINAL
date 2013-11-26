'''
ABOUT:
This program makes will make plots that help determine what is the best
pixfrac for the given scale.

DEPENDS:
multiprocessing
glob
pyfits
numpy
matplotlib
drizzlepac
stsci

AUTHOR:
Roberto J. Avila

HISTORY:
September 2013: Original script (Beta).

FUTURE IMPROVEMENTS:

USE:
python pixfrac_tester.py KERNEL PIXRANGE PIXRANGE
-nodriz option turns off drizzling and only makes statistics plot

You need a file in the running directory called <target>.reg
This file will tell the program what region of the weight image to use 
for statistics. The format of the file should be a single line that 
has the following information about the box limits:
x1 x2 y1 y2
'''

__author__ = 'Roberto J. Avila'
__version__ = 0.2


def drizzling(pixfrac):

    outname = fltr+'_'+kernel+'_'+str(pixfrac)
    logfile = outname

    #teal.unlearn('astrodrizzle')
    adriz.AstroDrizzle(imlist,\
                       output = outname,\
                       runfile = logfile,\
                       context = False,\
                       crbit = 4096,\
                       resetbits = 0,\
                       num_cores = 1,\
                       preserve = False,\
                       clean = True,\
                       in_memory = False,\
                       static = False,\
                       skysub = False,\
                       driz_separate = False,\
                       median = False,\
                       blot = False,\
                       driz_cr = False,\
                       driz_combine = True,\
                       final_kernel = kernel,\
                       final_pixfrac = pixfrac,\
                       final_fillval = 1e9,\
                       final_bits = 2400,\
                       final_wcs = True,\
                       final_rot = 0.0,\
                       final_outnx = None,\
                       final_outny = None,\
                       final_scale = scale)

    return

def wht_analysis():

    #Initializing variables
    whtlist = glob.glob('*'+kernel+'*_wht.fits')
    std_med = np.empty(len(whtlist),dtype=float)
    fraclist = np.empty(len(whtlist),dtype=float)

    #Loop that measures statistics, also some information gathering
    for i,im in enumerate(whtlist):

        print 'Gathering statistics for {}'.format(im)
        hdu = pyfits.open(im)
        hdr = hdu[0].header
        wht = hdu[0].data
        hdu.close()

        if i == 0:
            target = hdr['TARGNAME']
            scale = str(hdr['D001SCAL'])
            kern = hdr['D001KERN']
            if hdr['INSTRUME'] == 'ACS':
                nimg = hdr['NDRIZIM']/2
            else:
                nimg = hdr['NDRIZIM']

        #Read file that tells this function what region to measure the stats from
        x1,x2,y1,y2 = np.loadtxt(target+'.reg',usecols=(0,1,2,3),unpack=True,dtype=int)
        print x1,x2,y1,y2

        wht_std = np.std(wht[y1:y2,x1:x2])
        wht_med = np.median(wht[y1:y2,x1:x2])
        std_med[i] = wht_std/wht_med
        fraclist[i] = hdr['D001PIXF']
        
    #Plotting commands              
    plt.clf()
    plt.xlim(-0.025,1.025)
    plt.ylim(0.,1.05)
    plt.scatter(fraclist,std_med,s=50)
    plt.axhline(0.2,ls='--',lw=3,c='r')
    plt.xlabel('pixfrac',fontsize=18)
    plt.ylabel('rms/median',fontsize=18)
    plt.text(0.95,0.95,target,fontsize=16,horizontalalignment='right') 
    plt.text(0.95,0.90,fltr,fontsize=16,horizontalalignment='right')
    plt.text(0.95,0.85,'Scale='+str(scale)+'"',fontsize=16,horizontalalignment='right')
    plt.text(0.95,0.80,'Kernel='+kern,fontsize=16,horizontalalignment='right')
    plt.text(0.95,0.75,str(nimg)+' images',fontsize=16,horizontalalignment='right')
    plt.savefig(target+'_'+fltr+'_'+kern+'_'+str(scale)+'.png')

    return
        


if __name__ == "__main__":

    import multiprocessing
    import pyfits
    import glob
    import argparse
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from drizzlepac import astrodrizzle as adriz
    from stsci.tools import teal

    parser = argparse.ArgumentParser(description='Make drizzled images with\
                different pixfracs and measure statistics')
    parser.add_argument('kernel',type=str,default='square',\
                choices=['square','point','gaussian','turbo','tophat',\
                'lanczos3'],help='Kernel to use for drizzling.')
    parser.add_argument('pixrange',nargs=2,type=float,help=\
                        'Endpoints of  pixfracs range to test.')
    parser.add_argument('-nodriz','--no_driz',default=False,action='store_true',\
                        help='Turn off drizzling step. Used to only make plots')


    options = parser.parse_args()
    DODRIZ = not(options.no_driz)
    kernel=options.kernel
    pixfraclist = np.arange(options.pixrange[0],options.pixrange[1]+0.01,0.1)

    #Making an image list and grabbing information from the first in the list
    imlist = glob.glob('*fl?.fits')
    hdu = pyfits.open(imlist[0])
    hdr = hdu[0].header
    hdu.close()

    if hdr['INSTRUME'] == 'WFC3':
        fltr = hdr['FILTER']
        scale = 0.06
        imsize = 0.016 #size of image in GB

    if hdr['INSTRUME'] == 'ACS':
        tfltr = np.empty(2,dtype='|S5')
        tfltr[0] = hdr['FILTER1']
        tfltr[1] = hdr['FILTER2']
        fltr = tfltr[np.where(tfltr!='CLEAR')][0] 
        scale = 0.03
        imsize = 0.16 #size of image in GB

    #This is where the drizzling is done
    if DODRIZ:

        #Figure out memory usage
        #ACS file = 160MB, IR file = 16MB
        #8Kx8K drz image = 732MB x 2 (sci and wht images)
        #Do not exceed 32GB of memory usage
        mem_limit = 32.
        images_perproc = len(imlist)
        num_proc = len(pixfraclist)
        mem_usage = num_proc*(images_perproc*imsize+0.732*2)
        print 'Memory usage = {}GB'.format(mem_usage)

        if mem_usage <  mem_limit:

            print 'Using parallel processing'
            p = multiprocessing.Pool(len(pixfraclist))
            p.map(drizzling,pixfraclist)
            p.close()
            p.join()

        else:

            print 'Parallel processing turned off'
            for pixfraci in pixfraclist:

                drizzling(pixfraci)
        
        print '############################################################'
        print 'FINISHED DRIZZLING STEP'
        print '############################################################\n'

    
    #This is where the weight map analysis is done
    wht_analysis()

    print 'Finished'
