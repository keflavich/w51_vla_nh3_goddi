from __future__ import print_function
import os
import numpy as np
import pyspeckit
from astropy import units as u
from astropy import constants
from astropy import table
import radio_beam
from pyspeckit.spectrum.units import SpectroscopicAxis
from pyspeckit.spectrum.models import ammonia_constants,ammonia,ammonia_hf
import pylab as pl
amf = pyspeckit.spectrum.models.ammonia.ammonia_model_background()
from hf_only_model import hfonly_66_fixed_fitter, hfonly_fitter, sixsix_movinghf_fitter, sevenseven_movinghf_fitter

ammonia.freq_dict['ninenine'] = 27.47794e9
ammonia.freq_dict['tenten'] = 28.60475e9
ammonia.freq_dict['eleveneleven'] = 29.91449e9
ammonia.freq_dict['twelvetwelve'] = 31.42494e9
ammonia.freq_dict['thirteenthirteen'] = 33.15684e9

ammonia_constants.num_to_name[10] = 'ten'
ammonia_constants.num_to_name[11] = 'eleven'
ammonia_constants.num_to_name[12] = 'twelve'
ammonia_constants.num_to_name[13] = 'thirteen'

def load_spectrum(j=6, object='w51e2-tot',
                  npix=1,
                  headerfile='../W51-25GHzcont.map.image.fits',
                  fntemplate='spec-{object}-{j}{j2}-pb.txt',):
    # convert 6 -> 'six' for use below
    linename = ammonia_constants.num_to_name[j] * 2

    # the beam size is important for determining the brighntess conversion
    bm = radio_beam.Beam.from_fits_header(headerfile)

    # read in the file (dumped to text)
    xx,yy=np.loadtxt(fntemplate.format(object=object, j=j, j2=j if j < 10 else ""),
                     comments="#").T
    # we want MEAN flux, not SUM
    yy = yy / npix

    # x units are in km/s, reference is in Hz
    xarr = SpectroscopicAxis(xx*u.km/u.s,
                             refX=ammonia.freq_dict[linename]*u.Hz,
                             velocity_convention='radio')
    sp = pyspeckit.Spectrum(xarr=xarr, data=yy)
    sp.unit = 'Jy'
    # compute the typical noise over the -10 to 10 km/s region
    mean_error = sp.stats((-10,10))['std']
    sp.error[:] = mean_error
    sp.specname = '{0} {1}'.format(object, linename)

    # Copy the spectrum so we can convert it to Kelvins
    spK = sp.copy()
    jytok = ((1*u.Jy).to(u.K, u.brightness_temperature(bm, sp.xarr.refX)).value)
    sp.header['JYTOK'] = jytok
    spK.unit = 'K'
    spK.data = sp.data * jytok
    spK.error = sp.error * jytok
    mean_error_k = mean_error * jytok

    return sp, spK, mean_error, mean_error_k


def fit_highj_nh3_absorption(j=6, j2=6, object='w51e2-tot',
                             headerfile='../W51-25GHzcont.map.image.fits',
                             savepath='.',
                             fntemplate='spec-{object}-{j}{j2}-pb.txt',
                             background_tb=4.2e3):
    """
    Fit absorption lines of high-J NH3
    """

    linename = ammonia_constants.num_to_name[j] * 2
    
    sp, spK, mean_error, mean_error_k = load_spectrum(j=j, object=object, fntemplate=fntemplate)

    # create our ammonia fitter from the ammonia model that allows Tau as the
    # variable (instead of kinetic temperature)
    # the model comes in a few varieties: we've selected the one that assumes
    # there *is* a background, but we're fitting the continuum-subtracted
    # spectrum, so this one will accept the background as an input parameter
    amf = ammonia_hf.nh3_vtau[linename].background_contsub_fitter
    spK.specfit.Registry.add_fitter('ammonia_bg', amf, 5)

    sp.plotter()
    # Fit each component with a gaussian to determine the offsets between components
    sp.specfit(fittype='gaussian', guesses=[-0.15, 25, 3,
                                            -0.15, 32, 3,
                                            -0.30, 58, 6,
                                            -0.15, 84, 3,
                                            -0.15, 90, 3], quiet=True)
    # Store the fitted velocities and frequencies in convenient variables
    center = sp.specfit.parinfo[7].value
    voffs = u.Quantity([(((x-center)*u.km/u.s))
                        for x in sp.specfit.parinfo.values[1::3]])
    foffs = u.Quantity([((((x - center)*u.km/u.s) /
                          constants.c) * sp.xarr.refX).to(u.MHz)
                        for x in sp.specfit.parinfo.values[1::3]])
    print("Velocity Offsets: ",voffs)
    print("Frequency Offsets: ",foffs)
    models = (pyspeckit.spectrum.models.ammonia_constants.voff_lines_dict[linename][:3]
              +
              pyspeckit.spectrum.models.ammonia_constants.voff_lines_dict[linename][-2:])
    dfoffs = u.Quantity([(((x-center-m)*u.km/u.s))
                         for x,m in zip(sp.specfit.parinfo.values[1::3], models[::-1])])
    print("Difference between observed and modeled offsets: ",dfoffs)

    # plot the fits & residuals
    sp.specfit.plotresiduals(axis=sp.plotter.axis, clear=False,
                             yoffset=sp.data.min()*1.2, label=False)
    sp.plotter.axis.set_ylim(sp.data.min()*1.2 - 3*mean_error + sp.specfit.residuals.min(),
                             sp.data.max()*1.2 + 3*mean_error)

    sp.plotter.savefig(os.path.join(savepath,
                                    "{object}_nh3_{j}{j2}_eachcomponent.png".format(j=j,
                                                                                    j2=j2,
                                                                                    object=object)))

    T = True
    F = False
    #parnames=['Tbackground','Tex','tau','center','width'],
    # input guesses.  background_tb is fixed
    parinfo = amf.make_parinfo(params=[background_tb, 100, 2.5, 58, 0.5],
                               fixed=[T,F,F,F,F],)
    parinfo.TBACKGROUND0.limited = (False,False)
    parinfo.TEX0.mpmaxstep=1000
    print("DEBUG CHECK: Is parinfo set correctly? (should have background=fixed)\n",parinfo)

    # do the fits on the brightness temperature spectrum now (with ammonia model)
    # for the first fit, exclude the main line and fit only the hyperfines
    spK.plotter()
    spK.specfit.selectregion(exclude=[40,70], highlight=True)
    spK.specfit(fittype='ammonia_bg', parinfo=parinfo,
                #quiet=False, shh=False, verbose=True,
                veryverbose=False,
                reset_selection=False)

    spK.specfit.plot_fit()
    spK.specfit.plotresiduals(axis=spK.plotter.axis, clear=False,
                              yoffset=spK.data.min()*1.2, label=False)
    spK.plotter.axis.set_ylim(spK.data.min()*1.2 - 3*mean_error_k + spK.specfit.residuals.min(),
                              spK.data.max()*1.2 + 3*mean_error_k)
    spK.plotter.savefig("{object}_nh3_{j}{j2}_hyperfineonly.png".format(j=j,
                                                                        j2=j2,
                                                                        object=object))

    # second try: fit the main line and the hyperfines simultaneously
    spK.plotter()
    spK.specfit(fittype='ammonia_bg', parinfo=parinfo,
                #quiet=False, shh=False, verbose=True,
                veryverbose=False,
                reset_selection=True)
    spK.specfit.plotresiduals(axis=spK.plotter.axis, clear=False,
                              yoffset=spK.data.min()*1.2, label=False)
    spK.plotter.axis.set_ylim(spK.data.min()*1.2 - 3*mean_error_k + spK.specfit.residuals.min(),
                              spK.data.max()*1.2 + 3*mean_error_k)
    spK.plotter.savefig("{object}_nh3_{j}{j2}_fitwhole.png".format(j=j, j2=j2,
                                                                   object=object))

    # third try: all lines, but with two components
    parinfo = amf.make_parinfo(params=[background_tb, 200, 2.5, 58, 1.35,
                                       background_tb, 200, 2.5, 59,  0.5,],
                               fixed=[T,F,F,F,F,]*2,
                               npeaks=2)
    parinfo.TBACKGROUND0.limited = (False,False)
    parinfo.TBACKGROUND1.limited = (False,False)
    parinfo.TEX0.mpmaxstep=1000
    parinfo.TEX1.mpmaxstep=1000
    print("DEBUG CHECK: Is parinfo set correctly? (should have background=fixed and two sets of pars)\n",parinfo)
    spK.plotter()
    spK.specfit(fittype='ammonia_bg', parinfo=parinfo,
                #verbose=True, quiet=False, shh=False,
                veryverbose=False,
                reset_selection=True)
    spK.specfit.plotresiduals(axis=spK.plotter.axis, clear=False,
                              yoffset=spK.data.min()*1.2, label=False)
    spK.plotter.axis.set_ylim(spK.data.min()*1.2 - 3*mean_error_k + spK.specfit.residuals.min(),
                              spK.data.max()*1.2 + 3*mean_error_k)
    spK.plotter.axis.set_xlim(-30,200)
    spK.plotter.savefig("{object}_nh3_{j}{j2}_fitwhole_twocomp.png".format(j=j,
                                                                           j2=j2,
                                                                           object=object))

    pl.draw()
    pl.show()

    # TODO: refit with NH3 hyperfine model with velocity offsets and linewidths free
    # ??? optical depth reliable ???
    # examine Goddi 2015 on NGC 7538
    # -> column density (probably OK using simple LTE)
    # -> rotational temperature?

    spK.specfit.Registry.add_fitter('sixsix_movinghf', sixsix_movinghf_fitter(), 5)
    spK.specfit(fittype='sixsix_movinghf', guesses=[58, 30, 26.9, 31.2, 2])
    spK.specfit.plotresiduals(axis=spK.plotter.axis, clear=False,
                              yoffset=spK.data.min()*1.2, label=False)
    spK.plotter.axis.set_ylim(spK.data.min()*1.2 - 3*mean_error_k + spK.specfit.residuals.min(),
                              spK.data.max()*1.2 + 3*mean_error_k)
    spK.plotter.axis.set_xlim(-30,200)
    spK.plotter.savefig("{object}_nh3_{j}{j2}_fitwhole_movinghf.png".format(j=j,
                                                                            j2=j2,
                                                                            object=object))

    return sp, spK, (center, voffs, foffs, dfoffs)

def fit_all():
    emission_objects = ['w51e1', 'w51e2-em', ] # e1 is actually e8
    absorption_objects = ['w51e2-core', 'w51e2-south', 'w51e2-tot',]
    spectra = {"{object}_{j}{j2}".format(object=object, j=j, j2=j2):
               fit_highj_nh3_absorption(j=j, j2=j2, object=object)
               for object in absorption_objects
               for j,j2 in ((6,6),(7,7),)} #(9,9),(1,3),(1,0))}

    vals = [[key.split("_")[0],  key.split("_")[1], center]
            + list(voffs) + list(foffs) + list(dfoffs)
            for key,(sp,spK,(center, voffs, foffs, dfoffs))
            in spectra.items()]

    # confusing as heck but... trust me, k?
    columns = [table.Column(data=[d.value if hasattr(d,'value') else d for d in
                                  data],
                            name=name,
                            unit=data[0].unit if hasattr(data[0],'unit') else None)
               for data,name in zip(zip(*vals),
                                    ['Object', 'Line', 'Center',]+
                                    ['voffs{0}'.format(ii) for ii in range(5)]+
                                    ['foffs{0}'.format(ii) for ii in range(5)]+
                                    ['dfoffs{0}'.format(ii) for ii in range(5)],
                                   )]

    tbl = table.Table(columns)

    return spectra,tbl

if __name__ == "__main__":
    spectra,tbl = fit_all()
    print(tbl)

    fig = pl.figure(13)
    fig.clf()
    T=True
    F=False
    sp = spectra['w51e2-core_66'][1]
    mean_error = sp.data[-20:-2].std()
    sp.plotter.figure = fig
    sp.plotter.axis = pl.subplot(2,1,1)
    sp.plotter()
    sp.specfit.Registry.add_fitter('sixsix_movinghf', sixsix_movinghf_fitter(), 7)
    sp.specfit(fittype='sixsix_movinghf', guesses=[58, 83, 26.9, 31.2, 2, 1800, 8200], fixed=[F,F,F,F,F,T,T], annotate=False)
    sp.plotter.axis.set_xticklabels([])
    sp.specfit.plotresiduals(axis=sp.plotter.axis, clear=False,
                             yoffset=sp.data.min()*1.3, label=False)
    sp.plotter.axis.set_ylim(sp.data.min()*1.3 - 3*mean_error + sp.specfit.residuals.min(),
                             sp.data.max()*1.3 + 3*mean_error)
    sp.plotter.axis.set_xlim(0, 120)
    sp1 = spectra['w51e2-core_77'][1]
    mean_error = sp1.data[-20:-2].std()
    sp1.plotter.figure = fig
    sp1.plotter.axis = pl.subplot(2,1,2)
    sp1.plotter()
    sp1.specfit.Registry.add_fitter('sevenseven_movinghf', sevenseven_movinghf_fitter(), 7)
    sp1.specfit(fittype='sevenseven_movinghf', guesses=[58, 83, 27.3, 31.2, 1.32, 1800, 8200], fixed=[F,F,F,F,F,T,T], annotate=False)
    sp1.specfit.plotresiduals(axis=sp1.plotter.axis, clear=False,
                             yoffset=sp1.data.min()*1.3, label=False)
    sp1.plotter.axis.set_ylim(sp1.data.min()*1.3 - 3*mean_error + sp1.specfit.residuals.min(),
                             sp1.data.max()*1.3 + 3*mean_error)
    sp1.plotter.axis.set_xlim(0, 120)
    pl.subplots_adjust(hspace=0)
    fig.savefig("core_66_and_77_fits_freehf_fixedtex.png", bbox_inches='tight')
    fig.savefig("core_66_and_77_fits_freehf_fixedtex.eps", bbox_inches='tight')
