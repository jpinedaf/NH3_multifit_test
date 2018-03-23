import numpy as np
import matplotlib.pylab as plt
import matplotlib.pylab as pl

import pyspeckit
from pyspeckit.spectrum.models.ammonia import cold_ammonia_model
from multicube.subcube import SubCube
from astropy.io import fits
from extras import verify_parlims
import warnings

file_cube='./data/NGC1333_NH3_11_base_DR1.fits'
file_rms='./data/NGC1333_NH3_11_base_DR1_rms.fits'
file_tp='./data/NGC1333_NH3_11_base_DR1_Tpeak.fits'
file_mc_guess1c='./fits/NGC1333_NH3_MC_guess_npeaks1.fits'
file_mc_guess2c='./fits/NGC1333_NH3_MC_guess_npeaks2.fits'
file_best_fit_1c='./fits/NGC1333_NH3_MC_best_npeaks1.fits'
file_best_fit_2c='./fits/NGC1333_NH3_MC_best_npeaks2.fits'

download_data=False
do_1comp=False
do_2comp=False
do_inspect=True

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
snr_cut=15
multicore=40

if do_1comp:
    sc = SubCube(file_cube)
    sc.plot_spectrum(103,133)
    npars = 6 # for an ammonia model
    npeaks = 1
    #line_names = ['oneone', 'twotwo']
    line_names = ['oneone']
    fittype = 'cold_ammonia_x{}'.format(npeaks)
    fitmodel = cold_ammonia_model

    sc.specfit.Registry.add_fitter(fittype, npars=npars,
                                   function=fitmodel(line_names=line_names))
    sc.update_model(fittype)
    sc.specfit.fitter.npeaks = npeaks
    # parameters are Tk, Tex, log(N(NH3)), sigma_v, v_lsr
    minpars = [10., 3.0, 13.5, 0.05, 6.0, 0.5]
    maxpars = [10., 9.0, 15.0, 1.50, 9.5, 0.5]
    finesse = [1, 6, 6, 10, 10, 1]
    sc.make_guess_grid(minpars, maxpars, finesse)
    sc.generate_model(multicore=multicore)
    sc.snr_map=snr_map
    #sc.best_guess(sn_cut=snr_cut)

    import os
    if os.path.exists('fits') == False:
        os.mkdir('fits')

    # Do we search for the best guess or load from file?
    # (delete `file_mc_guess` to regenerate the models and guesses!)
    if os.path.isfile(file_mc_guess1c):
        sc.best_guesses = fits.getdata(file_mc_guess1c)
    else:
        sc.best_guess(sn_cut=snr_cut)
        fits.writeto(file_mc_guess1c, sc.best_guesses, hd)

    #sc.plot_spectrum(103,133)
    #sc.plotter.axis.plot(sc.xarr.value, sc.model_grid[sc._best_map[133,103]])
    #plt.imshow(sc.best_guesses[4,:,:], origin='lowest', cmap='RdYlBu_r')
    #plt.show()

    # fix Tkin and f2o for both l.o.s. components
    sc.fiteach_args['fixed'][0]=True
    sc.fiteach_args['fixed'][5]=True

    # make sure the parameter limits are all within the bounds
    verify_parlims(sc.best_guesses, sc.fiteach_args,
                   npeaks=npeaks, npars=npars)

    sc.fiteach(fittype = fittype,
        guesses   = sc.best_guesses,
        multicore = multicore,
        errmap    = rms,
        verbose   = 0,
        **sc.fiteach_args)
    sc.write_fit(file_best_fit_1c,overwrite=True)

if do_2comp:
    sc2 = SubCube(file_cube)
    npeaks = 2 # number of l.o.s. components
    npars = 6 # for an ammonia model
    #line_names = ['oneone', 'twotwo']
    line_names = ['oneone']
    fittype = 'cold_ammonia_x{}'.format(npeaks)
    fitmodel = cold_ammonia_model

    sc2.specfit.Registry.add_fitter(fittype, npars=npars,
                                    function=fitmodel(line_names=line_names))
    sc2.update_model(fittype)
    sc2.specfit.fitter.npeaks = npeaks

    # parameters are Tk, Tex, log(N(NH3)), sigma_v, v_lsr
    minpars2= [10., 3.0, 13.5, 0.05, 6.0, 0.5]+[10., 3.0, 13.5, 0.05, 6.0, 0.5]
    maxpars2= [10., 9.0, 15.0, 1.50, 9.5, 0.5]+[10., 9.0, 15.0, 1.50, 9.5, 0.5]
    finesse2= [1, 1, 4, 6, 8, 1] + [1, 1, 4, 6, 8, 1]
    sc2.make_guess_grid(minpars2, maxpars2, finesse2)
    sc2.generate_model(multicore=multicore)
    sc2.snr_map=snr_map
    #sc2.get_snr_map()

    import os
    if os.path.exists('fits') == False:
        os.mkdir('fits')

    # Do we search for the best guess or load from file?
    # (delete `file_mc_guess` to regenerate the models and guesses!)
    if os.path.isfile(file_mc_guess2c):
        sc2.best_guesses = fits.getdata(file_mc_guess2c)
    else:
        sc2.best_guess(sn_cut=snr_cut)
        fits.writeto(file_mc_guess2c, sc2.best_guesses, hd)

    # for the highest S/N pixel on the map
    j, i = np.unravel_index(np.nanargmax(sc2.snr_map), sc2.snr_map.shape)
    #i=106 # 130
    #j=159 # 144
    sc2.plot_spectrum(i,j)
    ij_model = sc2.specfit.get_full_model(pars=sc2.best_guesses[:, j, i])
    sc2.plotter.axis.plot(sc2.xarr.value, ij_model)
    plt.show()

    # fix Tkin and f2o for both l.o.s. components
    sc2.fiteach_args['fixed'][0]=True
    sc2.fiteach_args['fixed'][5]=True
    sc2.fiteach_args['fixed'][6]=True
    sc2.fiteach_args['fixed'][11]=True

    # make sure the parameter limits are all within the bounds
    verify_parlims(sc2.best_guesses, sc2.fiteach_args,
                   npeaks=npeaks, npars=npars)

    sc2.fiteach(fittype   = sc2.fittype,
        guesses   = sc2.best_guesses,
        multicore = multicore,
        errmap    = rms,
        verbose   = 0,
        **sc2.fiteach_args)

    # vsokolov: not sure why this is needed, for some reason parinfo does not
    # have a correct npeaks after fiteach with npeaks > 1...
    parinfo, _ = sc2.specfit.fitter._make_parinfo(npeaks=npeaks)
    sc2.specfit.parinfo = sc2.specfit.fitter.parinfo = parinfo
    sc2.write_fit(file_best_fit_2c,overwrite=True)

if do_inspect:
    i=109
    j=141

    # Load cube 1
    sc = SubCube(file_cube)
    npars = 6 # for an ammonia model
    npeaks = 1
    line_names = ['oneone']
    fittype = 'cold_ammonia_x{}'.format(npeaks)
    fitmodel = cold_ammonia_model
    sc.specfit.Registry.add_fitter(fittype, npars=npars,
                                   function=fitmodel(line_names=line_names))
    sc.update_model(fittype)
    sc.specfit.fitter.npeaks = npeaks
    # Load Data 2
    sc2 = SubCube(file_cube)
    npeaks = 2 # number of l.o.s. components
    fittype = 'cold_ammonia_x{}'.format(npeaks)
    fitmodel = cold_ammonia_model
    sc2.specfit.Registry.add_fitter(fittype, npars=npars,
                                    function=fitmodel(line_names=line_names))
    sc2.update_model(fittype)
    sc2.specfit.fitter.npeaks = npeaks
    #
    # plotting
    #
    best_fit1c= fits.getdata(file_best_fit_1c)
    best_fit2c= fits.getdata(file_best_fit_2c)
    sc2.plot_spectrum(i,j, axis=pl.subplot(2,1,1))
    ij_model1 = sc.specfit.get_full_model(pars=best_fit1c[0:6, j, i])
    ij_model2 = sc.specfit.get_full_model(pars=best_fit2c[0:6, j, i])
    ij_model3 = sc.specfit.get_full_model(pars=best_fit2c[6:12, j, i])
    ij_model4 = sc2.specfit.get_full_model(pars=best_fit2c[0:12, j, i])
    sc2.plotter.axis.plot(sc2.xarr.value, ij_model1, color='r')
    sc2.plotter.axis.plot(sc2.xarr.value, ij_model2, color='g')
    sc2.plotter.axis.plot(sc2.xarr.value, ij_model3, color='b')
    sc2.plotter.axis.plot(sc2.xarr.value, ij_model4, color='orange')

    pl.subplot(2,1,2).plot(sc2.xarr.value, sc2.data-ij_model1, color='r')
    #pl.subplot(2,1,2).plot(sc2.xarr.value, sc2.data-ij_model2, color='g')
    #pl.subplot(2,1,2).plot(sc2.xarr.value, sc2.data-ij_model3, color='b')
    pl.subplot(2,1,2).plot(sc2.xarr.value, sc2.data-ij_model4, color='orange')
