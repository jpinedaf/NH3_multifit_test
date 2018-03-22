import numpy as np
import matplotlib.pylab as plt

import pyspeckit
from multicube.subcube import SubCube
from astropy.io import fits
import warnings

file_cube='./data/NGC1333_NH3_11_base_DR1.fits'
file_rms='./data/NGC1333_NH3_11_base_DR1_rms.fits'
file_tp='./data/NGC1333_NH3_11_base_DR1_Tpeak.fits'
file_mc_guess='./fits/NGC1333_NH3_MC_guess.fits'

download_data=False
do_1comp=False
do_2comp=True

if download_data:
    from astropy.utils.data import download_file
    import os, shutil
    if os.path.exists('data') == False:
        os.mkdir('data')
    # get rms file
    if os.path.isfile(file_rms) == False:
        tmp_rms=download_file('https://dataverse.harvard.edu/api/access/datafile/3008213')
        shutil.move(tmp_rms, file_rms)
    # get cube
    if os.path.isfile(file_cube) == False:
        tmp_cube=download_file('https://dataverse.harvard.edu/api/access/datafile/2902348')
        shutil.move(tmp_cube, file_cube)
    # create Peak temperature file
    if os.path.isfile(file_tp) == False:
        hd=fits.getheader(file_rms)
        hd['BUNIT']='K'
        cube=fits.getdata(file_cube)
        tp_data=cube.max(axis=0)
        fits.writeto(file_tp, tp_data, hd)

hd=fits.getheader(file_rms)
tp= fits.getdata(file_tp)
rms= fits.getdata(file_rms)
snr_map=tp/rms
snr_cut=4

if do_1comp:
    sc = SubCube(file_cube)
    sc.plot_spectrum(103,133)
    sc.update_model('cold_ammonia')
    # parameters are Tk, Tex, log(N(NH3)), sigma_v, v_lsr
    minpars = [10., 3.0, 13.5, 0.05, 6.0, 0.5]
    maxpars = [10., 9.0, 15.0, 1.50, 9.5, 0.5]
    finesse = [1, 6, 6, 10, 10, 1]
    sc.make_guess_grid(minpars, maxpars, finesse)
    sc.generate_model()
    #sc.get_snr_map()
    sc.snr_map=snr_map
    sc.best_guess(sn_cut=snr_cut)
    #sc.plot_spectrum(103,133)
    #sc.plotter.axis.plot(sc.xarr.value, sc.model_grid[sc._best_map[133,103]])
    #plt.imshow(sc.best_guesses[4,:,:], origin='lowest', cmap='RdYlBu_r')
    #plt.show()

multicore=40
if do_2comp:
    # parameters are Tk, Tex, log(N(NH3)), sigma_v, v_lsr
    sc2= SubCube(file_cube)
    npeaks = 2 # number of l.o.s. components
    line_names = ['oneone', 'twotwo']
    line_names = ['oneone']
    fittype_fmt = 'cold_ammonia_x{}'
    from pyspeckit.spectrum.models.ammonia import cold_ammonia_model
    npars = 6 # for an ammonia model
    fitmodel = cold_ammonia_model

    
    sc2.specfit.Registry.add_fitter(fittype_fmt.format(npeaks), npars=npars,
                                    function=fitmodel(line_names=line_names))
    sc2.update_model(fittype_fmt.format(npeaks))
    sc2.specfit.fitter.npeaks = npeaks
    #
    minpars2= [10., 3.0, 13.5, 0.05, 6.0, 0.5]+[10., 3.0, 13.5, 0.05, 6.0, 0.5]
    maxpars2= [10., 9.0, 15.0, 1.50, 9.5, 0.5]+[10., 9.0, 15.0, 1.50, 9.5, 0.5]
    finesse2= [1, 1, 4, 8,10, 1] + [1, 1, 4, 8,10, 1]
    sc2.make_guess_grid(minpars2, maxpars2, finesse2)
    sc2.generate_model(multicore=multicore)
    sc2.snr_map=snr_map
    #sc2.get_snr_map()
    sc2.best_guess(sn_cut=snr_cut)
    
    import os
    if os.path.exists('fits') == False:
        os.mkdir('fits')
    #
    # Do we search for the best guess or load from file?
    if os.path.isfile(file_mc_guess):
        sc2._best_map = fits.getdata( file_mc_guess)
    else:
        fits.writeto( file_mc_guess, sc2._best_map, hd)

    #
    i=106 # 130
    j=159 # 144
    sc2.plot_spectrum(i,j)
    sc2.plotter.axis.plot(sc2.xarr.value, sc2.model_grid[sc2._best_map[j,i]])
    plt.show()
    #
    from multicube.astro_toolbox import get_ncores

    sc2.fiteach_args['fixed'][0]=True
    sc2.fiteach_args['fixed'][5]=True
    sc2.fiteach_args['fixed'][6]=True
    sc2.fiteach_args['fixed'][11]=True

    sc2.fiteach(fittype   = sc2.fittype,
        guesses   = sc2.best_guesses, 
        multicore = multicore,
        errmap    = rms,
        verbose   = 0,
        **sc2.fiteach_args)
    sc2.write_fit('fit_cube_filename.fits',overwrite=True)

