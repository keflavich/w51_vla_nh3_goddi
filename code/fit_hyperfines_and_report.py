"""
Code for populating table 3

Run this in a directory that contains the relevant spectrum text files
"""
import os
import numpy as np
import pyspeckit
import goddi_nh3_fits
from hf_only_model import hfonly_fitter, sixsix_movinghf_fitter
from astropy import units as u
from astropy import constants
from astroquery.splatalogue import Splatalogue
import pylab as pl

# we're going to make a lot of figures.  Close some.
for ii in pl.get_fignums():
    pl.close(ii)

savepath = '/Users/adam/Dropbox/w51_goddi_share/spec-casa/figures/'
resultpath = '/Users/adam/Dropbox/w51_goddi_share/fit_results.txt'

outf = open(resultpath, 'w')

distance = 5.41*u.kpc
XNH3 = 1e-7

areas = {'e2e': np.pi*0.45*u.arcsec*0.4459*u.arcsec,
         'e2abs': np.pi*0.45*u.arcsec*0.40*u.arcsec,
         'e8': np.pi*2.1707*u.arcsec*1.6168*u.arcsec,
         'e2nw': np.pi*0.450*u.arcsec*0.400*u.arcsec,
         'e2clump': (10.87-0.58)*u.arcsec**2,
        }

freqs = {'66': 25.05596*u.GHz,
         '77': 25.71514*u.GHz,
         '99': 27.47794*u.GHz,
         '1010': 28.60474*u.GHz,
         '1313': 33.15685*u.GHz,
        }
eupper = {'66': 409,
          '77': 539,
          '99': 853,
          '1010': 1037,
          '1313': 1693,
         }
degeneracy = {'66': 2*(2*6+1),
              '77': 1*(2*7+1),
              '99': 2*(2*9+1),
              '1010': 1*(2*10+1),
              '1313': 1*(2*13+1),
             }

nh3tables = {line:Splatalogue.query_lines(freqs[line]*0.99999,
                                          freqs[line]*1.00001,
                                          line_strengths=['ls1','ls2','ls3','ls4','ls5'],
                                          chemical_name=' NH3 ')
             for line in freqs}
linestrengths = {'66':float(nh3tables['66'][nh3tables['66']['Linelist']=='JPL']['S<sub>ij</sub>&#956;<sup>2</sup> (D<sup>2</sup>)'])*1e-18*u.esu*u.cm,
                 '77':float(nh3tables['77'][nh3tables['77']['Linelist']=='JPL']['S<sub>ij</sub>&#956;<sup>2</sup> (D<sup>2</sup>)'])*1e-18*u.esu*u.cm,
                 '99':float(nh3tables['99'][nh3tables['99']['Linelist']=='JPL']['S<sub>ij</sub>&#956;<sup>2</sup> (D<sup>2</sup>)'])*1e-18*u.esu*u.cm,
                 '1010':float(nh3tables['1010'][nh3tables['1010']['Linelist']=='JPL']['S<sub>ij</sub>&#956;<sup>2</sup> (D<sup>2</sup>)'])*1e-18*u.esu*u.cm,
                 '1313':float(nh3tables['1313'][nh3tables['1313']['Linelist']=='JPL']['S<sub>ij</sub>&#956;<sup>2</sup> (D<sup>2</sup>)'])*1e-18*u.esu*u.cm,
                }

FWHM = np.sqrt(8*np.log(2))

def tbl_vals(pi):
    """
    Compute the average amplitude, average width, and respective errors from the hyperfine fits
    """
    amp = (pi[2] + pi[5])/2.
    amperr = (pi[2].error + pi[5].error)/2.
    sigma = (pi[3] + pi[6])/2.
    sigmaerr = (pi[3].error + pi[6].error)/2.
    integral = amp * (2*np.pi)**0.5 *(sigma)
    int_error = (((amperr/amp)**2 + (sigmaerr/sigma)**2)*integral**2)**0.5

    return amp, amperr, sigma*FWHM, sigmaerr*FWHM, integral, int_error, pi[0].value, pi[0].error

def integral_hf(pi):
    return tbl_vals(pi)[4:6]

def tbl_vals_gaussian(pi):
    """
    Compute the average amplitude, average width, and respective errors from the hyperfine fits
    """
    amp = pi[0].value
    amperr = pi[0].error
    sigma = pi[2].value
    sigmaerr = pi[2].error
    integral = amp * (2*np.pi)**0.5 *(sigma)
    int_error = (((amperr/amp)**2 + (sigmaerr/sigma)**2)*integral**2)**0.5

    return amp, amperr, sigma*FWHM, sigmaerr*FWHM, integral, int_error, pi[1].value, pi[1].error

def integral_ml(pi):
    return tbl_vals_gaussian(pi)[4:6]

F = False
T = True

sp6e8,spK6e8,_,_ = goddi_nh3_fits.load_spectrum(6, object='w51e8', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp6e8.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp7e8,spK7e8,_,_ = goddi_nh3_fits.load_spectrum(7, object='w51e8', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp7e8.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp9e8,spK9e8,_,_ = goddi_nh3_fits.load_spectrum(9, object='w51e8', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp9e8.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp10e8,spK10e8,_,_ = goddi_nh3_fits.load_spectrum(10, object='w51e8', headerfile='/Users/adam/work/w51/goddi/W51-29GHzcont.map.image.fits')
sp10e8.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp13e8,spK13e8,_,_ = goddi_nh3_fits.load_spectrum(13, object='w51e8', headerfile='/Users/adam/work/w51/goddi/W51-33GHzcont.map.image.fits')
sp13e8.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)

sp6e8.plotter()
sp7e8.plotter()
sp9e8.plotter()
sp10e8.plotter()
sp13e8.plotter()
limits = [(50, 65), (25,30), (0, 1), (0.1, 6), (30, 35), (0, 1), (0.1, 6)]
limited = [(True,True)]*7
tied = ['']*7
tied[6] = 'p[3]'
tied[5] = 'p[2]'
sp6e8.specfit(fittype='hfonly', guesses=[58, 26.9, 0.02, 2.0, 31.4, 0.02, 2.0], fixed=[F,T,F,F,T,F,F], limited=limited, limits=limits, tied=tied)
sp7e8.specfit(fittype='hfonly', guesses=[58, 27.3, 0.02, 2.0, 31.2, 0.02, 2.0], fixed=[F,T,F,F,T,F,F], limited=limited, limits=limits, tied=tied)
sp9e8.specfit(fittype='hfonly', guesses=[58, 27.0, 0.02, 2.0, 30.1, 0.02, 2.0], fixed=[F,T,F,F,T,F,F], limited=limited, limits=limits, tied=tied)
sp6e8.specfit.plot_fit(show_hyperfine_components=True)
sp7e8.specfit.plot_fit(show_hyperfine_components=True)
sp9e8.specfit.plot_fit(show_hyperfine_components=True)
hfpi6e8 = sp6e8.specfit.parinfo
hfpi7e8 = sp7e8.specfit.parinfo
hfpi9e8 = sp9e8.specfit.parinfo
print("e8", file=outf)
print("hf 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp6e8.specfit.parinfo))), file=outf)
print("hf 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp7e8.specfit.parinfo))), file=outf)
"""
6-6 peak=0.0552312979281 +/- 0.00339213398418, fwhm=12.0092504285 +/- 0.346442710066, integral=3.60073813128 +/- 0.265490554173, center=59.8513485763 +/- 0.165718109489
7-7 peak=0.0343580596483 +/- 0.00211133669018, fwhm=7.08292886544 +/- 0.477361456933, integral=0.779163995174 +/- 0.0883610449159, center=59.9800528242 +/- 0.29151363486
"""
sp6e8.plotter.savefig(os.path.join(savepath, 'sp6e8_hf.png'))
sp7e8.plotter.savefig(os.path.join(savepath, 'sp7e8_hf.png'))
sp9e8.plotter.savefig(os.path.join(savepath, 'sp9e8_hf.png'))
sp10e8.plotter.savefig(os.path.join(savepath, 'sp10e8_hf.png'))
sp13e8.plotter.savefig(os.path.join(savepath, 'sp13e8_hf.png'))


sp6e8.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp7e8.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp9e8.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp10e8.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp13e8.specfit(fittype='gaussian', guesses=[1, 58, 1])
mlpi6e8 = sp6e8.specfit.parinfo
mlpi7e8 = sp7e8.specfit.parinfo
mlpi9e8 = sp9e8.specfit.parinfo
mlpi10e8 = sp9e8.specfit.parinfo
mlpi13e8 = sp9e8.specfit.parinfo
print("gaussian 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp6e8.specfit.parinfo))), file=outf)
print("gaussian 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp7e8.specfit.parinfo))), file=outf)
print("gaussian 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp9e8.specfit.parinfo))), file=outf)
print("gaussian 10-10 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp10e8.specfit.parinfo))), file=outf)
print("gaussian 13-13 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp13e8.specfit.parinfo))), file=outf)
"""
gaussian 6-6 peak=0.278613851186 +/- 0.00324112286397, fwhm=12.3888928639 +/- 0.166414956455, integral=19.3304597881 +/- 0.430595034151, center=59.0587207674 +/- 0.0706699251669
gaussian 7-7 peak=0.19046626119 +/- 0.00375662817502, fwhm=11.1985631899 +/- 0.255043061619, integral=10.7973493377 +/- 0.407787818476, center=59.7258654571 +/- 0.108306812353
gaussian 9-9 peak=0.0829717146021 +/- 0.00468072809892, fwhm=11.1981635835 +/- 0.729471923465, integral=4.70325110054 +/- 0.508070309565, center=59.6499937385 +/- 0.309778206795
"""



# e2w: absorption
sp6e2w,spK6e2w,_,_ = goddi_nh3_fits.load_spectrum(6, object='w51e2w', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp6e2w.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp7e2w,spK7e2w,_,_ = goddi_nh3_fits.load_spectrum(7, object='w51e2w', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp7e2w.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp9e2w,spK9e2w,_,_ = goddi_nh3_fits.load_spectrum(9, object='w51e2w', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp9e2w.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp10e2w,spK10e2w,_,_ = goddi_nh3_fits.load_spectrum(10, object='w51e2w', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp10e2w.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp13e2w,spK13e2w,_,_ = goddi_nh3_fits.load_spectrum(13, object='w51e2w', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp13e2w.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)

sp6e2w.plotter()
sp7e2w.plotter()
sp9e2w.plotter()
limits = [(50, 65), (25,30), (-10, 0), (0.1, 6), (30, 35), (-10, 0), (0.1, 6)]
limited = [(True,True)]*7
tied = ['']*7
tied[6] = 'p[3]'
tied[5] = 'p[2]'
sp6e2w.specfit(fittype='hfonly', guesses=[58, 26.9, -0.1, 2.0, 31.4, -0.1, 2.0], fixed=[F,F,F,F,F,F,F], limited=limited, limits=limits)
sp7e2w.specfit(fittype='hfonly', guesses=[58, 27.3, -0.1, 2.0, 31.2, -0.1, 2.0], fixed=[F,F,F,F,F,F,F], limited=limited, limits=limits)
sp9e2w.specfit(fittype='hfonly', guesses=[58, 27.0, -0.1, 2.0, 30.1, -0.1, 2.0], fixed=[F,T,F,F,T,F,F], limited=limited, limits=limits, tied=tied)
sp6e2w.specfit.plot_fit(show_hyperfine_components=True)
sp7e2w.specfit.plot_fit(show_hyperfine_components=True)
sp9e2w.specfit.plot_fit(show_hyperfine_components=True)
hfpi6e2w = sp6e2w.specfit.parinfo
hfpi7e2w = sp7e2w.specfit.parinfo
hfpi9e2w = sp9e2w.specfit.parinfo
sp6e2w.plotter.savefig(os.path.join(savepath,'sp6e2w_hf.png'))
sp7e2w.plotter.savefig(os.path.join(savepath,'sp7e2w_hf.png'))
sp9e2w.plotter.savefig(os.path.join(savepath,'sp9e2w_hf.png'))

print("e2w", file=outf)
print("hf 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp6e2w.specfit.parinfo))), file=outf)
print("hf 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp7e2w.specfit.parinfo))), file=outf)
print("hf 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp9e2w.specfit.parinfo))), file=outf)
sp6e2w.specfit(fittype='gaussian', guesses=[-1, 58, 1])
sp7e2w.specfit(fittype='gaussian', guesses=[-1, 58, 1])
sp9e2w.specfit(fittype='gaussian', guesses=[-1, 58, 1])
sp10e2w.specfit(fittype='gaussian', guesses=[-1, 58, 1])
sp13e2w.specfit(fittype='gaussian', guesses=[-1, 58, 1])
mlpi6e2w = sp6e2w.specfit.parinfo
mlpi7e2w = sp7e2w.specfit.parinfo
mlpi9e2w = sp9e2w.specfit.parinfo
mlpi10e2w = sp9e2w.specfit.parinfo
mlpi13e2w = sp9e2w.specfit.parinfo
print("gaussian 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp6e2w.specfit.parinfo))), file=outf)
print("gaussian 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp7e2w.specfit.parinfo))), file=outf)
print("gaussian 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp9e2w.specfit.parinfo))), file=outf)
print("gaussian 10-10 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp10e2w.specfit.parinfo))), file=outf)
print("gaussian 13-13 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp13e2w.specfit.parinfo))), file=outf)


# e2e aka e2east
sp6e2e,spK6e2e,_,_ = goddi_nh3_fits.load_spectrum(6, object='w51e2e', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp6e2e.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
spK6e2e.specfit.Registry.add_fitter('sixsix',sixsix_movinghf_fitter(),7)
sp7e2e,spK7e2e,_,_ = goddi_nh3_fits.load_spectrum(7, object='w51e2e', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp7e2e.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp9e2e,spK9e2e,_,_ = goddi_nh3_fits.load_spectrum(9, object='w51e2e', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp9e2e.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp10e2e,spK10e2e,_,_ = goddi_nh3_fits.load_spectrum(10, object='w51e2e', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp10e2e.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp13e2e,spK13e2e,_,_ = goddi_nh3_fits.load_spectrum(13, object='w51e2e', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp13e2e.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)

sp6e2e.plotter()
spK6e2e.plotter()
sp7e2e.plotter()
sp9e2e.plotter()
# sp6e2e might be optically thick (in hyperfines, in emission!)
sp6e2e.specfit(fittype='hfonly', guesses=[58, 26.9, 0.02, 2.0, 31.4, 0.02, 2.0], fixed=[F,T,F,F,T,F,F], tied=tied)
spK6e2e.specfit(fittype='sixsix', guesses=[58, 100, 26.9, 31.4, 2.0, 1200, 0], fixed=[F,F,T,T,F,T,T])
spK6e2e.specfit.plot_fit(show_hyperfine_components=True)
sp7e2e.specfit(fittype='hfonly', guesses=[58, 27.3, 0.02, 2.0, 31.2, 0.02, 2.0], fixed=[F,T,F,F,T,F,F], tied=tied)
sp9e2e.specfit(fittype='hfonly', guesses=[58, 27.0, 0.02, 2.0, 30.1, 0.02, 2.0], fixed=[F,T,F,F,T,F,F], tied=tied)
sp6e2e.specfit.plot_fit(show_hyperfine_components=True)
sp7e2e.specfit.plot_fit(show_hyperfine_components=True)
sp9e2e.specfit.plot_fit(show_hyperfine_components=True)
hfpi6e2e = sp6e2e.specfit.parinfo
hfpi7e2e = sp7e2e.specfit.parinfo
hfpi9e2e = sp9e2e.specfit.parinfo
sp6e2e.plotter.savefig(os.path.join(savepath,'sp6e2e_hf.png'))
spK6e2e.plotter.savefig(os.path.join(savepath,'spK6e2e_hf.png'))
sp7e2e.plotter.savefig(os.path.join(savepath,'sp7e2e_hf.png'))
sp9e2e.plotter.savefig(os.path.join(savepath,'sp9e2e_hf.png'))

print("e2e", file=outf)
print("hf 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp6e2e.specfit.parinfo))), file=outf)
print("hf 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp7e2e.specfit.parinfo))), file=outf)
print("hf 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp9e2e.specfit.parinfo))), file=outf)
sp6e2e.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp7e2e.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp9e2e.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp10e2e.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp13e2e.specfit(fittype='gaussian', guesses=[1, 58, 1])
mlpi6e2e = sp6e2e.specfit.parinfo
mlpi7e2e = sp7e2e.specfit.parinfo
mlpi9e2e = sp9e2e.specfit.parinfo
mlpi10e2e = sp9e2e.specfit.parinfo
mlpi13e2e = sp9e2e.specfit.parinfo
print("gaussian 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp6e2e.specfit.parinfo))), file=outf)
print("gaussian 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp7e2e.specfit.parinfo))), file=outf)
print("gaussian 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp9e2e.specfit.parinfo))), file=outf)
print("gaussian 10-10 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp10e2e.specfit.parinfo))), file=outf)
print("gaussian 13-13 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp13e2e.specfit.parinfo))), file=outf)


# e2nw, e2northwest
sp6e2nw,spK6e2nw,_,_ = goddi_nh3_fits.load_spectrum(6, object='w51e2nw', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp6e2nw.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp7e2nw,spK7e2nw,_,_ = goddi_nh3_fits.load_spectrum(7, object='w51e2nw', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp7e2nw.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp9e2nw,spK9e2nw,_,_ = goddi_nh3_fits.load_spectrum(9, object='w51e2nw', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp9e2nw.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp10e2nw,spK10e2nw,_,_ = goddi_nh3_fits.load_spectrum(10, object='w51e2nw', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp10e2nw.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp13e2nw,spK13e2nw,_,_ = goddi_nh3_fits.load_spectrum(13, object='w51e2nw', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp13e2nw.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)

sp6e2nw.plotter()
sp7e2nw.plotter()
sp9e2nw.plotter()
sp6e2nw.specfit(fittype='hfonly', guesses=[58, 26.9, 0.02, 2.0, 31.4, 0.02, 2.0], fixed=[F,T,F,F,T,F,F], tied=tied)
sp7e2nw.specfit(fittype='hfonly', guesses=[58, 27.3, 0.02, 2.0, 31.2, 0.02, 2.0], fixed=[F,T,F,F,T,F,F], tied=tied)
limits = [(50, 65), (25,30), (0, 1), (0.1, 6), (30, 35), (0, 1), (0.1, 6)]
sp9e2nw.specfit(fittype='hfonly', guesses=[58, 27.0, 0.02, 2.0, 30.1, 0.02, 2.0], fixed=[F,T,F,F,T,F,F], tied=tied, limited=limited, limits=limits)
sp6e2nw.specfit.plot_fit(show_hyperfine_components=True)
sp7e2nw.specfit.plot_fit(show_hyperfine_components=True)
sp9e2nw.specfit.plot_fit(show_hyperfine_components=True)
hfpi6e2nw = sp6e2nw.specfit.parinfo
hfpi7e2nw = sp7e2nw.specfit.parinfo
hfpi9e2nw = sp9e2nw.specfit.parinfo
sp6e2nw.plotter.savefig(os.path.join(savepath,'sp6e2nw_hf.png'))
sp7e2nw.plotter.savefig(os.path.join(savepath,'sp7e2nw_hf.png'))
sp9e2nw.plotter.savefig(os.path.join(savepath,'sp9e2nw_hf.png'))

print("e2nw", file=outf)
print("hf 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp6e2nw.specfit.parinfo))), file=outf)
print("hf 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp7e2nw.specfit.parinfo))), file=outf)
print("hf 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp9e2nw.specfit.parinfo))), file=outf)
sp6e2nw.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp7e2nw.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp9e2nw.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp10e2nw.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp13e2nw.specfit(fittype='gaussian', guesses=[1, 58, 1])
mlpi6e2nw = sp6e2nw.specfit.parinfo
mlpi7e2nw = sp7e2nw.specfit.parinfo
mlpi9e2nw = sp9e2nw.specfit.parinfo
mlpi10e2nw = sp9e2nw.specfit.parinfo
mlpi13e2nw = sp9e2nw.specfit.parinfo
print("gaussian 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp6e2nw.specfit.parinfo))), file=outf)
print("gaussian 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp7e2nw.specfit.parinfo))), file=outf)
print("gaussian 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp9e2nw.specfit.parinfo))), file=outf)
print("gaussian 10-10 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp10e2nw.specfit.parinfo))), file=outf)
print("gaussian 13-13 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp13e2nw.specfit.parinfo))), file=outf)


# e2clump
sp6e2clump,spK6e2clump,_,_ = goddi_nh3_fits.load_spectrum(6, object='w51e2-clump', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp6e2clump.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp7e2clump,spK7e2clump,_,_ = goddi_nh3_fits.load_spectrum(7, object='w51e2-clump', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp7e2clump.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp9e2clump,spK9e2clump,_,_ = goddi_nh3_fits.load_spectrum(9, object='w51e2-clump', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp9e2clump.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp10e2clump,spK10e2clump,_,_ = goddi_nh3_fits.load_spectrum(10, object='w51e2-clump', headerfile='/Users/adam/work/w51/goddi/W51-29GHzcont.map.image.fits')
sp10e2clump.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp13e2clump,spK13e2clump,_,_ = goddi_nh3_fits.load_spectrum(13, object='w51e2-clump', headerfile='/Users/adam/work/w51/goddi/W51-33GHzcont.map.image.fits')
sp13e2clump.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)

sp6e2clump.plotter()
sp7e2clump.plotter()
sp9e2clump.plotter()
limits = [(50, 65), (25,30), (0, 1), (0.1, 6), (30, 35), (0, 1), (0.1, 6)]
limited = [(True,True)]*7
sp6e2clump.specfit(fittype='hfonly', guesses=[58, 26.9, 0.02, 2.0, 31.4, 0.02, 2.0], fixed=[F,T,F,F,T,F,F], limited=limited, limits=limits, tied=tied)
sp7e2clump.specfit(fittype='hfonly', guesses=[58, 27.3, 0.02, 2.0, 31.2, 0.02, 2.0], fixed=[F,T,F,F,T,F,F], limited=limited, limits=limits, tied=tied)
sp9e2clump.specfit(fittype='hfonly', guesses=[58, 27.0, 0.02, 2.0, 30.1, 0.02, 2.0], fixed=[F,T,F,F,T,F,F], limited=limited, limits=limits, tied=tied)
hfpi6e2clump = sp6e2clump.specfit.parinfo
hfpi7e2clump = sp7e2clump.specfit.parinfo
hfpi9e2clump = sp9e2clump.specfit.parinfo
sp6e2clump.specfit.plot_fit(show_hyperfine_components=True)
sp7e2clump.specfit.plot_fit(show_hyperfine_components=True)
sp9e2clump.specfit.plot_fit(show_hyperfine_components=True)
sp6e2clump.plotter.savefig(os.path.join(savepath,'sp6e2clump_hf.png'))
sp7e2clump.plotter.savefig(os.path.join(savepath,'sp7e2clump_hf.png'))
sp9e2clump.plotter.savefig(os.path.join(savepath,'sp9e2clump_hf.png'))

print("e2clump", file=outf)
print("hf 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp6e2clump.specfit.parinfo))), file=outf)
print("hf 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp7e2clump.specfit.parinfo))), file=outf)
print("hf 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp9e2clump.specfit.parinfo))), file=outf)
sp6e2clump.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp7e2clump.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp9e2clump.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp10e2clump.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp13e2clump.specfit(fittype='gaussian', guesses=[1, 58, 1])
mlpi6e2clump = sp6e2clump.specfit.parinfo
mlpi7e2clump = sp7e2clump.specfit.parinfo
mlpi9e2clump = sp9e2clump.specfit.parinfo
mlpi10e2clump = sp9e2clump.specfit.parinfo
mlpi13e2clump = sp9e2clump.specfit.parinfo
print("gaussian 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp6e2clump.specfit.parinfo))), file=outf)
print("gaussian 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp7e2clump.specfit.parinfo))), file=outf)
print("gaussian 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp9e2clump.specfit.parinfo))), file=outf)
print("gaussian 10-10 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp10e2clump.specfit.parinfo))), file=outf)
print("gaussian 13-13 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp13e2clump.specfit.parinfo))), file=outf)


outf.close()



def full_integral(pi_hf, pi_ml):
    """
    Combine the hyperfine and main line integral
    """
    if pi_hf is not None:
        i1,ei1 = integral_hf(pi_hf)
    else:
        i1,ei1 = 0,0

    i2,ei2 = integral_ml(pi_ml)

    return (i1+i2), (ei1**2+ei2**2)**0.5

def Nu_thin_hightex(flux, line_strength, freq, fillingfactor=1.0, tau=None):
    """
    Optically-thin-ish approximation for the column density of the upper state
    of a given line assuming T_ex >> T_bg and T_ex >> h nu
    """
    assert flux.unit.is_equivalent(u.K*u.km/u.s)
    assert line_strength.unit.is_equivalent(u.esu*u.cm)
    k = constants.k_B
    term1 = (3*k/(8*np.pi**2 * freq * line_strength**2))
    term5 = flux.to(u.K*u.km/u.s) / fillingfactor
    term6 = 1 if tau is None else tau/(1-np.exp(-tau))
    return (term1*term5*term6).to(u.cm**-2)

flux_e8 = {'66': (full_integral(hfpi6e8, mlpi6e8)*u.Jy).to(u.K, u.brightness_temperature(areas['e8'], freqs['66']))*u.km/u.s,
           '77': (full_integral(hfpi7e8, mlpi7e8)*u.Jy).to(u.K, u.brightness_temperature(areas['e8'], freqs['77']))*u.km/u.s,
           '99': (full_integral(hfpi9e8, mlpi9e8)*u.Jy).to(u.K, u.brightness_temperature(areas['e8'], freqs['99']))*u.km/u.s,
           '1010': (full_integral(None, mlpi10e8)*u.Jy).to(u.K, u.brightness_temperature(areas['e8'], freqs['1010']))*u.km/u.s,
           '1313': (full_integral(None, mlpi13e8)*u.Jy).to(u.K, u.brightness_temperature(areas['e8'], freqs['1313']))*u.km/u.s,
          }
flux_e2e = {'66': (full_integral(hfpi6e2e, mlpi6e2e)*u.Jy).to(u.K, u.brightness_temperature(areas['e2e'], freqs['66']))*u.km/u.s,
            '77': (full_integral(hfpi7e2e, mlpi7e2e)*u.Jy).to(u.K, u.brightness_temperature(areas['e2e'], freqs['77']))*u.km/u.s,
            '99': (full_integral(hfpi9e2e, mlpi9e2e)*u.Jy).to(u.K, u.brightness_temperature(areas['e2e'], freqs['99']))*u.km/u.s,
            '1010': (full_integral(None, mlpi10e2e)*u.Jy).to(u.K, u.brightness_temperature(areas['e2e'], freqs['1010']))*u.km/u.s,
            '1313': (full_integral(None, mlpi13e2e)*u.Jy).to(u.K, u.brightness_temperature(areas['e2e'], freqs['1313']))*u.km/u.s,
          }
taue2e = {'66': 100,
            '77': 77,
            '99': 43,
            '1010': None,
            '1313': None,
           }
taue8 = {'66':27,
         '77':33,
         '99':None,
         '1010': None,
         '1313': None,
           }
# if you want to exclude 1010 from the fit, uncomment this
#freqs.pop('1010')

Nu_e8 = {line: Nu_thin_hightex(flux_e8[line], linestrengths[line], freqs[line]) for line in freqs}
Nu_e2e = {line: Nu_thin_hightex(flux_e2e[line], linestrengths[line], freqs[line]) for line in freqs}
Nu_e8_taucorr = {line: Nu_thin_hightex(flux_e8[line], linestrengths[line], freqs[line], tau=taue8[line]) for line in freqs}
Nu_e2e_taucorr = {line: Nu_thin_hightex(flux_e2e[line], linestrengths[line], freqs[line], tau=taue2e[line]) for line in freqs}


flux_e2nw = {'66': (full_integral(hfpi6e2nw, mlpi6e2nw)*u.Jy).to(u.K, u.brightness_temperature(areas['e2nw'], freqs['66']))*u.km/u.s,
            '77': (full_integral(hfpi7e2nw, mlpi7e2nw)*u.Jy).to(u.K, u.brightness_temperature(areas['e2nw'], freqs['77']))*u.km/u.s,
            '99': (full_integral(hfpi9e2nw, mlpi9e2nw)*u.Jy).to(u.K, u.brightness_temperature(areas['e2nw'], freqs['99']))*u.km/u.s,
#            '1010': (full_integral(None, mlpi10e2nw)*u.Jy).to(u.K, u.brightness_temperature(areas['e2nw'], freqs['1010']))*u.km/u.s,
#            '1313': (full_integral(None, mlpi13e2nw)*u.Jy).to(u.K, u.brightness_temperature(areas['e2nw'], freqs['1313']))*u.km/u.s,
          }
taue2nw = {'66': None,
            '77': None,
            '99': None,
            '1010': None,
            '1313': None,
           }
Nu_e2nw = {line: Nu_thin_hightex(flux_e2nw[line], linestrengths[line], freqs[line]) for line in freqs if line in flux_e2nw}


flux_e2clump = {'66': (full_integral(hfpi6e2clump, mlpi6e2clump)*u.Jy).to(u.K, u.brightness_temperature(areas['e2clump'], freqs['66']))*u.km/u.s,
            '77': (full_integral(hfpi7e2clump, mlpi7e2clump)*u.Jy).to(u.K, u.brightness_temperature(areas['e2clump'], freqs['77']))*u.km/u.s,
            '99': (full_integral(hfpi9e2clump, mlpi9e2clump)*u.Jy).to(u.K, u.brightness_temperature(areas['e2clump'], freqs['99']))*u.km/u.s,
#            '1010': (full_integral(None, mlpi10e2clump)*u.Jy).to(u.K, u.brightness_temperature(areas['e2clump'], freqs['1010']))*u.km/u.s,
#            '1313': (full_integral(None, mlpi13e2clump)*u.Jy).to(u.K, u.brightness_temperature(areas['e2clump'], freqs['1313']))*u.km/u.s,
          }
taue2clump = {'66': None,
            '77': None,
            '99': None,
            '1010': None,
            '1313': None,
           }
Nu_e2clump = {line: Nu_thin_hightex(flux_e2clump[line], linestrengths[line], freqs[line]) for line in freqs if line in flux_e2clump}

# http://spec.jpl.nasa.gov/ftp/pub/catalog/doc/d017002.pdf
B0=(298117.0*u.MHz)
C0=(186726.0*u.MHz)

from astropy import modeling

def fit_tex(eupper, nupperoverg):
    """
    Fit the Boltzmann diagram
    """
    model = modeling.models.Linear1D()
    fitter = modeling.fitting.LevMarLSQFitter()
    result = fitter(model, eupper, np.log(nupperoverg))
    tex = -1./result.slope*u.K

    #calculating partition function Q_rot; use eqn 15.48 tools of radio astronomy
    # = sqrt( (pi (kT)^3) / h^3 B^2 C )
    Q_rot= np.sqrt(np.pi * (constants.k_B*tex)**3 / (constants.h**3 * B0**2 * C0)).decompose()

    Ntot = np.exp(result.intercept + np.log(Q_rot)) * u.cm**-2

    print(("Tex={0}, Ntot={1}, Q_rot={2}".format(tex, Ntot, Q_rot)))

    return Ntot, tex, result.slope, result.intercept

import pylab as pl
pl.matplotlib.rc_file('pubfiguresrc')

eups = np.array([eupper[line] for line in freqs])

fig = pl.figure(1)
fig.clf()
ax = pl.gca()
print("e8")
a,b,c = ax.errorbar(x=[eupper[line] for line in freqs],
            y=[Nu_e8[line][0].value*degeneracy[line] for line in freqs],
            yerr=[Nu_e8[line][1].value*degeneracy[line] for line in freqs],
            marker='s',
            linestyle='none')
Ntot,tex,slope,intcpt = fit_tex([eupper[line] for line in freqs], [Nu_e8[line][0].value*degeneracy[line] for line in freqs])
ax.plot(eups, np.exp(eups*slope + intcpt), color=a.get_color(), label='$T_{{ex}}={0:0.1f}$\n$N(\\mathrm{{NH}}_3)={1:0.1e}$ cm$^{{-2}}$'.format(tex, Ntot.value))
print(("Mass lower limit: {0}".format((Ntot/XNH3 * areas['e8']*distance**2 * 2.8*u.Da).to(u.M_sun, u.dimensionless_angles()))))
a,b,c = ax.errorbar(x=[eupper[line] for line in freqs],
            y=[Nu_e8_taucorr[line][0].value*degeneracy[line] for line in freqs],
            yerr=[Nu_e8_taucorr[line][1].value*degeneracy[line] for line in freqs],
            marker='s',
            linestyle='none')
Ntot,tex,slope,intcpt = fit_tex([eupper[line] for line in freqs],
                                [Nu_e8_taucorr[line][0].value*degeneracy[line] for line in freqs])
ax.plot(eups, np.exp(eups*slope + intcpt), color=a.get_color(), label='$T_{{ex}}={0:0.1f}$\n$N(\\mathrm{{NH}}_3)={1:0.1e}$ cm$^{{-2}}$'.format(tex, Ntot.value))
ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_ylabel("$N_u / g$ [cm$^{-2}$]")
ax.set_xlabel("$E_u$ [K]")
ax.set_title('e8')
pl.legend(loc='best')
pl.savefig(os.path.join(savepath, 'e8_rotational_diagram_fit.png'))

print("e2e")
fig = pl.figure(2)
fig.clf()
ax = pl.gca()
a,b,c = ax.errorbar(x=[eupper[line] for line in freqs],
            y=[Nu_e2e[line][0].value*degeneracy[line] for line in freqs],
            yerr=[Nu_e2e[line][1].value*degeneracy[line] for line in freqs],
            marker='s',
            linestyle='none')
Ntot,tex,slope,intcpt = fit_tex([eupper[line] for line in freqs], [Nu_e2e[line][0].value*degeneracy[line] for line in freqs])
ax.plot(eups, np.exp(eups*slope + intcpt), color=a.get_color(), label='$T_{{ex}}={0:0.1f}$\n$N(\\mathrm{{NH}}_3)={1:0.1e}$ cm$^{{-2}}$'.format(tex, Ntot.value))
print(("Mass lower limit: {0}".format((Ntot/XNH3 * areas['e2e']*distance**2 * 2.8*u.Da).to(u.M_sun, u.dimensionless_angles()))))
a,b,c = ax.errorbar(x=[eupper[line] for line in freqs],
            y=[Nu_e2e_taucorr[line][0].value*degeneracy[line] for line in freqs],
            yerr=[Nu_e2e_taucorr[line][1].value*degeneracy[line] for line in freqs],
            marker='s',
            linestyle='none')
Ntot,tex,slope,intcpt = fit_tex([eupper[line] for line in freqs],
                                [Nu_e2e_taucorr[line][0].value*degeneracy[line] for line in freqs])
ax.plot(eups, np.exp(eups*slope + intcpt), color=a.get_color(), label='$T_{{ex}}={0:0.1f}$\n$N(\\mathrm{{NH}}_3)={1:0.1e}$ cm$^{{-2}}$'.format(tex, Ntot.value))
ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_ylabel("$N_u / g$ [cm$^{-2}$]")
ax.set_xlabel("$E_u$ [K]")
ax.set_title('e2e')
pl.legend(loc='best')
pl.savefig(os.path.join(savepath, 'e2e_rotational_diagram_fit.png'))

print("e2nw")
fig = pl.figure(3)
fig.clf()
ax = pl.gca()
a,b,c = ax.errorbar(x=[eupper[line] for line in freqs if line in Nu_e2nw],
            y=[Nu_e2nw[line][0].value*degeneracy[line] for line in freqs if line in Nu_e2nw],
            yerr=[Nu_e2nw[line][1].value*degeneracy[line] for line in freqs if line in Nu_e2nw],
            marker='s',
            linestyle='none')
Ntot,tex,slope,intcpt = fit_tex([eupper[line] for line in freqs if line in Nu_e2nw], [Nu_e2nw[line][0].value*degeneracy[line] for line in freqs if line in Nu_e2nw])
ax.plot(eups, np.exp(eups*slope + intcpt), color=a.get_color(), label='$T_{{ex}}={0:0.1f}$\n$N(\\mathrm{{NH}}_3)={1:0.1e}$ cm$^{{-2}}$'.format(tex, Ntot.value))
print(("Mass lower limit: {0}".format((Ntot/XNH3 * areas['e2nw']*distance**2 * 2.8*u.Da).to(u.M_sun, u.dimensionless_angles()))))
ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_ylabel("$N_u / g$ [cm$^{-2}$]")
ax.set_xlabel("$E_u$ [K]")
ax.set_title('e2nw')
pl.legend(loc='best')
pl.savefig(os.path.join(savepath, 'e2nw_rotational_diagram_fit.png'))


print("e2clump")
fig = pl.figure(4)
fig.clf()
ax = pl.gca()
a,b,c = ax.errorbar(x=[eupper[line] for line in freqs if line in Nu_e2clump],
            y=[Nu_e2clump[line][0].value*degeneracy[line] for line in freqs if line in Nu_e2clump],
            yerr=[Nu_e2clump[line][1].value*degeneracy[line] for line in freqs if line in Nu_e2clump],
            marker='s',
            linestyle='none')
Ntot,tex,slope,intcpt = fit_tex([eupper[line] for line in freqs if line in Nu_e2clump], [Nu_e2clump[line][0].value*degeneracy[line] for line in freqs if line in Nu_e2clump])
ax.plot(eups, np.exp(eups*slope + intcpt), color=a.get_color(), label='$T_{{ex}}={0:0.1f}$\n$N(\\mathrm{{NH}}_3)={1:0.1e}$ cm$^{{-2}}$'.format(tex, Ntot.value))
print(("Mass lower limit: {0}".format((Ntot/XNH3 * areas['e2clump']*distance**2 * 2.8*u.Da).to(u.M_sun, u.dimensionless_angles()))))
ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_ylabel("$N_u / g$ [cm$^{-2}$]")
ax.set_xlabel("$E_u$ [K]")
ax.set_title('e2clump')
pl.legend(loc='best')
pl.savefig(os.path.join(savepath, 'e2clump_rotational_diagram_fit.png'))

pl.draw()
pl.show()


