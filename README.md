WW-pipeline-FINAL
=================
Latest scripts and configuration files that will be used to align and drizzle images. The order for calling the scripts is as follows:

1. *run_make_crclean*
    * _Purpose_: construct crclean images for flts located in asn file (default), or any specified flts.
    * _Output_: [root]_crclean.fits; crclean.log (default)
    * _Options_:  
               [-a asn filename] [-log  AD log filename] [-c  AD # cores]  
               [-na "no asn" switch to use all flts in cwd] [-im  flt images to use when -na switch is on]

2. *run_sex_crclean*
    * _Purpose_: construct SExtractor catalogs for crclean images using config [instrum].sex.crclean.config.
    * makes trimmed catalog by keeping the brightest objects (100 is default) with reliable SEx flags (0 or 2).
    * creates ds9 region files for full and trimmed catalog (via the "cat2reg" gawk script).
    * creates external ascii file "catfile.sex" that will be used in step 3.
    * _Output_: *.sex.all = full catalog; *.sex = trimmed catalog; *.reg = ds9 region file (two);  
               wfc3.sex.xml = SE logfile; catfile.sex = tweakreg input file (step 3)
    * _Options_: [-im images] [-max max objects to keep in trimmed catalog] [-th SE threshold]
    * _WARNING_: expects default.param, default.conv, default.nnw, cat2reg, and  
               [instrum].sex.crclean.config in cwd.

    -- OR --

   *run_find_crclean* (temporary name holder)
    * _Purpose_: construct ImageFind catalogs for each crlean image.
    * removes sources located within 50 pixels of image edge.
    * creates external ascii file "catfile.coo" that will be used in step 3.
    * Output: *.coo = trimmed ImageFind catalog; catfile.coo = tweakreg input

3. *run_tweakreg_flt*
    * _Purpose_: align flt images with tweakreg using external catalogs (e.g., SExtractor,ImageFind, hst2align).
    * _Output_: tweakreg.log (default name); imlist.dat = flt image list; *.match = catalog of matched objects
    * _Options_:  
               [-im images] [-cf catfile name] [-log tweakreg logfile name] [-rim ref image] [-rcat refcat]  
               [-xc catalog xcol] [-yc catalog ycol] [-rxc refcat xcol] [-ryc refcat ycol]

4. *run_make_tweakfigs*
    * _Purpose_: construct diagnostic diagrams (delta-x/y and vectorgram) using tweakreg *.match files.
    * _Output_: [root]_xyresids.pdf & [root]_xyvector.pdf
    * _Options_: [-cat tweakreg *fit.match catalogs]
    * _WARNING_: expects corresponding flt images in same directory.

5. *run_drizzle*
    * _Purpose_: create drizzled image from a stack of aligned flt images.
    * _Output_: drizzled sci, wht, and ctx images (root=[filter]_[whttype]); astrodrizzle.log (default name);  
               "imlist.dat" - flt image list
    * _Options_:  
         [-im images] [-c # cores] [-pf AD final_pixfrac] [-ps AD final_scale]  
         [-fb AD final_bits] [-fk AD final_kernel] [-wht AD final_wht_type]  
         [-ra AD final_ra] [-dec AD final_dec] [-nx AD final_outnx] [-ny AD final_outny]
    * _WARNING_: recommend running it twice, using ivm and exp weighting. This gives ivm and exposure map that can be used for weighting and object flags, respectively, in the next step.

6. *run_sex_drz*
    * _Purpose_: construct SExtractor catalogs for drizzled images using config named [instrum].sex.drz.config.
    * makes trimmed catalog same as step2 above (i.e., select on flags & flux), but also optionally selects objects using minimum CLASS_STAR and maximum FWHM_IMAGE values.
    * creates ds9 region files for full and trimmed catalog (via the "cat2reg" gawk script).
    * _Output_: *.sex.all = full catalog; *.sex = trimmed catalog; *.reg = ds9 region file (2);  
               wfc3.sex.xml = SE logfile; catfile.sex = tweakreg input file (step 7);  
               [filter]_rms.fits (optional;RMS map); [filter]_expflag.fits (optional; exposure flag map)
    * _Options_:  
               [-im images] [-dim detection image] [-ivm IVM weight image] [-exp EXP weight image] [-th SE threshold]  
               [-max max objects to keep in trimmed catalog] [-cs min CLASS_STAR for trimmed catalog]  
               [-f max FWHM_IMAGE for trimmed catalog]
    * _WARNING_:  
      (a) requires default.param, default.conv, default.nnw, cat2reg, & [instrum].sex.drz.config in cwd.  
      (b) will create a flag map based on the exposure image (if provided) to exclude bad sources from trimmed catalog.

7. *run_tweakreg_drz*
    * _Purpose_: align drz to drz images (using script similar to version for aligning flt images).
    * _Output_: tweakreg.log; imlist.dat = list of flt images to align; *fit.match = catalog of matched objects
    * _Options_:  
               [-im images] [-cf catfile name] [-log tweakreg logfile name] [-rim ref image] [-rcat ref catalog]  
               [-xc catalog x column] [-yc catalog y column] [-rxc refcat x column] [-ryc refcat y column]
    * _WARNING_: main must include if/else statement:  
         IF   filter == reference, then refim/refcat correspond to ground-based data,  
         ELSE refim/refcat correspond to our chosen reference filter.

8. *run_make_tweakfigs*
    * _Purpose_: constructs tweakreg diagnostic diagrams using drz to drz *fit.match files (same script used in step 4).
    * _Output_: [drzroot]_xyresids.pdf & [drzroot]_xyvector.pdf
    * _Options_: [-cat tweakreg *fit.match catalogs]
    * _WARNING_: expects corresponding drz images in cwd.

9. *run_tweakback*
    * _Purpose_: apply wcs solution from aligned drizzled image to constituent flt images.
    * _Output_: edits WCS in flt images. No new files are created.
    * _Options_: [-drz drizzled image]
    * _WARNING_: expects that the constituent flt/flc images are located in cwd.


10. *pixfrac_tester*
    * _Purpose_: Makes drizzled images with different pixfracs and then plots statistics to select "best pixfrac."
    * _Output_: drizzled sci and wht images for each pixfrac; A png file with the rms/median statistics.
    * _Options_: [-nodriz] Turns off drizzling and only makes stats plot
    * _WARNING_: expects corresponding fl? images in cwd. Also needs a file called \<targetname\>.reg which is used to define the box region where stats will be measured


11. *run_drizzle*
    * _Purpose_: drizzle tweakback'ed flt images to a common survey footprint using "best pixfrac" determined in previous step (same script as step 5).
    * _Output_: drizzled sci, wht, and ctx images (root=[filter]_[whttype]); astrodrizzle.log (default name); "imlist.dat" - flt image list
    * _Options_:  
         [-im images] [-c # cores] [-pf AD final_pixfrac] [-ps AD final_scale]  
         [-fb AD final_bits] [-fk AD final_kernel] [-wht AD final_wht_type]  
         [-ra AD final_ra] [-dec AD final_dec] [-nx AD final_outnx] [-ny AD final_outny]
    * _WARNING_: we recommend running it twice using ivm and exp weighting - allows for rms and flag maps as input to SExtractor.

12. *run_sex_drz*
    * _Purpose_: construct SExtractor catalogs for drizzled footprint images using config named [instrum].sex.drz.config.
    * makes trimmed catalog same as step2 above (i.e., select on flags & flux), but also optionally selects objects using minimum CLASS_STAR and maximum FWHM_IMAGE values.
    * creates ds9 region files for full and trimmed catalog (via the "cat2reg" gawk script).
    * _Output_: *.sex.all = full catalog; *.sex = trimmed catalog; *.reg = ds9 region file (2);
               wfc3.sex.xml = SE logfile; catfile.sex = tweakreg input file (step 7);
               [filter]_rms.fits (optional;RMS map); [filter]_expflag.fits (optional; exposure flag map)
    * _Options_:  
               [-im images] [-dim detection image] [-ivm IVM weight image] [-exp EXP weight image] [-th SE threshold]
               [-max max objects to keep in trimmed catalog] [-cs min CLASS_STAR for trimmed catalog]
               [-f max FWHM_IMAGE for trimmed catalog]
    * _WARNING_:  
      (a) requires default.param, default.conv, default.nnw, cat2reg, & [instrum].sex.drz.config in cwd.

---

A typical sequence of calls from the shell assuming that the fl? images we are aligning/drizzling are located in cwd.
The parameter choices seen below assume we are aligning at least ~10 IR flts - adjustments will be necessary depending on instrument, filter, and the number of flt/flc images. Here, we also ignore bookkeeping to save original
copies and organize output files, etc.

    >python run_make_crclean.py  
    >python run_sex_crclean.py  
    >python run_tweakreg_flt.py  
    >python run_make_tweakfigs.py  
    >python run_drizzle.py -ps 0.06             # IR (for few flts use native -pf 1.0 -ps 0.13")  
    >python run_drizzle.py -ps 0.06 -wht EXP  
    >python run_sex_drz.py -im [filt]_ivm_drz_sci.fits -ivm [filt]_ivm_drz_wht.fits -exp [filt]_exp_drz_wht.fits  
    >python run_tweakreg_drz.py -im [filt]_ivm_drz_sci.fits -rim [reffilt]_ivm_drz_sci.fits -rcat [reffilt]_ivm_drz_sci.sex  
    >python run_make_tweakfigs.py  
    >python run_tweakback.py -drz [filt]_ivm_drz_sci.fits  
    >python pixfrac_tester.py 0.3 1.0  
    >python run_drizzle.py -pf best_pixfrac -ps 0.06 -ra 196.0 -dec 27.0 -nx 8000 -ny 8000
    >python run_drizzle.py -pf best_pixfrac -ps 0.06 -ra 196.0 -dec 27.0 -nx 8000 -ny 8000 -wht EXP  
    >python run_sex_drz.py -im [filt]_ivm_drz_sci.fits -dim [reffilt]_ivm_drz_sci.fits -ivm [filt]_ivm_drz_wht.fits  
    -exp [filt]_exp_drz_wht.fits -th 3

