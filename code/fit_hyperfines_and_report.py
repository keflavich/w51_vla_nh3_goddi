"""
Code for populating table 3
"""
import numpy as np
import pyspeckit
import goddi_nh3_fits
from hf_only_model import hfonly_fitter
from astropy import units as u
from astropy import constants
from astroquery.splatalogue import Splatalogue

distance = 5.41*u.kpc
XNH3 = 1e-7

areas = {'e2emi': np.pi*0.45*u.arcsec*0.4459*u.arcsec,
         'e2abs': np.pi*0.45*u.arcsec*0.40*u.arcsec,
         'e8': np.pi*2.1707*u.arcsec*1.6168*u.arcsec,
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

sp6e2emi,spK6e2emi,_,_ = goddi_nh3_fits.load_spectrum(6, object='w51e2-proto', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp6e2emi.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp7e2emi,spK7e2emi,_,_ = goddi_nh3_fits.load_spectrum(7, object='w51e2-proto', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp7e2emi.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp9e2emi,spK9e2emi,_,_ = goddi_nh3_fits.load_spectrum(9, object='w51e2-proto', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp9e2emi.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp10e2emi,spK10e2emi,_,_ = goddi_nh3_fits.load_spectrum(10, object='w51e2-proto', headerfile='/Users/adam/work/w51/goddi/W51-29GHzcont.map.image.fits')
sp10e2emi.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp13e2emi,spK13e2emi,_,_ = goddi_nh3_fits.load_spectrum(13, object='w51e2-proto', headerfile='/Users/adam/work/w51/goddi/W51-33GHzcont.map.image.fits')
sp13e2emi.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)

F = False
T = True
sp6e2emi.specfit(fittype='hfonly', guesses=[58, 26.9, 0.02, 2.0, 31.4, 0.02, 2.0], fixed=[F,F,F,F,F,F,F])
sp7e2emi.specfit(fittype='hfonly', guesses=[58, 26.9, 0.02, 2.0, 31.4, 0.02, 2.0], fixed=[F,F,F,F,F,F,F])
sp9e2emi.specfit(fittype='hfonly', guesses=[58, 27.0, 0.02, 2.0, 30.1, 0.02, 2.0], fixed=[F,T,F,F,T,F,F])
hfpi6e2emi = sp6e2emi.specfit.parinfo
hfpi7e2emi = sp7e2emi.specfit.parinfo
hfpi9e2emi = sp9e2emi.specfit.parinfo

print "e2emi"
print "hf 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp6e2emi.specfit.parinfo)))
print "hf 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp7e2emi.specfit.parinfo)))
print "hf 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp9e2emi.specfit.parinfo)))
sp6e2emi.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp7e2emi.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp9e2emi.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp10e2emi.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp13e2emi.specfit(fittype='gaussian', guesses=[1, 58, 1])
mlpi6e2emi = sp6e2emi.specfit.parinfo
mlpi7e2emi = sp7e2emi.specfit.parinfo
mlpi9e2emi = sp9e2emi.specfit.parinfo
mlpi10e2emi = sp9e2emi.specfit.parinfo
mlpi13e2emi = sp9e2emi.specfit.parinfo
print "gaussian 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp6e2emi.specfit.parinfo)))
print "gaussian 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp7e2emi.specfit.parinfo)))
print "gaussian 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp9e2emi.specfit.parinfo)))
print "gaussian 10-10 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp10e2emi.specfit.parinfo)))
print "gaussian 13-13 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp13e2emi.specfit.parinfo)))

sp6e8,spK6e8,_,_ = goddi_nh3_fits.load_spectrum(6, object='w51e8', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp6e8.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp7e8,spK7e8,_,_ = goddi_nh3_fits.load_spectrum(7, object='w51e8', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp7e8.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp9e8,spK9e8,_,_ = goddi_nh3_fits.load_spectrum(9, object='w51e8', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp9e8.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp10e8,spK10e8,_,_ = goddi_nh3_fits.load_spectrum(10, object='w51e2-proto', headerfile='/Users/adam/work/w51/goddi/W51-29GHzcont.map.image.fits')
sp10e8.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp13e8,spK13e8,_,_ = goddi_nh3_fits.load_spectrum(13, object='w51e2-proto', headerfile='/Users/adam/work/w51/goddi/W51-33GHzcont.map.image.fits')
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
sp6e8.specfit(fittype='hfonly', guesses=[58, 26.9, 0.02, 2.0, 31.4, 0.02, 2.0], fixed=[F,F,F,F,F,F,F], limited=limited, limits=limits)
sp7e8.specfit(fittype='hfonly', guesses=[58, 26.9, 0.02, 2.0, 31.4, 0.02, 2.0], fixed=[F,F,F,F,F,F,F], limited=limited, limits=limits, tied=tied)
sp9e8.specfit(fittype='hfonly', guesses=[58, 27.0, 0.02, 2.0, 30.1, 0.02, 2.0], fixed=[F,T,F,F,T,F,F], limited=limited, limits=limits, tied=tied)
hfpi6e8 = sp6e8.specfit.parinfo
hfpi7e8 = sp7e8.specfit.parinfo
hfpi9e8 = sp9e8.specfit.parinfo
print "e8"
print "hf 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp6e8.specfit.parinfo)))
print "hf 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp7e8.specfit.parinfo)))
"""
6-6 peak=0.0552312979281 +/- 0.00339213398418, fwhm=12.0092504285 +/- 0.346442710066, integral=3.60073813128 +/- 0.265490554173, center=59.8513485763 +/- 0.165718109489
7-7 peak=0.0343580596483 +/- 0.00211133669018, fwhm=7.08292886544 +/- 0.477361456933, integral=0.779163995174 +/- 0.0883610449159, center=59.9800528242 +/- 0.29151363486
"""
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
print "gaussian 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp6e8.specfit.parinfo)))
print "gaussian 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp7e8.specfit.parinfo)))
print "gaussian 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp9e8.specfit.parinfo)))
print "gaussian 10-10 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp10e8.specfit.parinfo)))
print "gaussian 13-13 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp13e8.specfit.parinfo)))
"""
gaussian 6-6 peak=0.278613851186 +/- 0.00324112286397, fwhm=12.3888928639 +/- 0.166414956455, integral=19.3304597881 +/- 0.430595034151, center=59.0587207674 +/- 0.0706699251669
gaussian 7-7 peak=0.19046626119 +/- 0.00375662817502, fwhm=11.1985631899 +/- 0.255043061619, integral=10.7973493377 +/- 0.407787818476, center=59.7258654571 +/- 0.108306812353
gaussian 9-9 peak=0.0829717146021 +/- 0.00468072809892, fwhm=11.1981635835 +/- 0.729471923465, integral=4.70325110054 +/- 0.508070309565, center=59.6499937385 +/- 0.309778206795
"""


sp6e2abs,spK6e2abs,_,_ = goddi_nh3_fits.load_spectrum(6, object='w51e2-tot', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp6e2abs.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp7e2abs,spK7e2abs,_,_ = goddi_nh3_fits.load_spectrum(7, object='w51e2-tot', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp7e2abs.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp9e2abs,spK9e2abs,_,_ = goddi_nh3_fits.load_spectrum(9, object='w51e2-tot', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp9e2abs.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp10e2abs,spK10e2abs,_,_ = goddi_nh3_fits.load_spectrum(10, object='w51e2-tot', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp10e2abs.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp13e2abs,spK13e2abs,_,_ = goddi_nh3_fits.load_spectrum(13, object='w51e2-tot', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp13e2abs.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)

sp6e2abs.plotter()
sp7e2abs.plotter()
sp9e2abs.plotter()
limits = [(50, 65), (25,30), (-10, 0), (0.1, 6), (30, 35), (-10, 0), (0.1, 6)]
limited = [(True,True)]*7
tied = ['']*7
tied[6] = 'p[3]'
tied[5] = 'p[2]'
sp6e2abs.specfit(fittype='hfonly', guesses=[58, 26.9, -0.1, 2.0, 31.4, -0.1, 2.0], fixed=[F,F,F,F,F,F,F], limited=limited, limits=limits)
sp7e2abs.specfit(fittype='hfonly', guesses=[58, 26.9, -0.1, 2.0, 31.4, -0.1, 2.0], fixed=[F,F,F,F,F,F,F], limited=limited, limits=limits)
sp9e2abs.specfit(fittype='hfonly', guesses=[58, 27.0, -0.1, 2.0, 30.1, -0.1, 2.0], fixed=[F,T,F,F,T,F,F], limited=limited, limits=limits, tied=tied)
hfpi6e2abs = sp6e2abs.specfit.parinfo
hfpi7e2abs = sp7e2abs.specfit.parinfo
hfpi9e2abs = sp9e2abs.specfit.parinfo
print "e2abs"
print "hf 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp6e2abs.specfit.parinfo)))
print "hf 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp7e2abs.specfit.parinfo)))
print "hf 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp9e2abs.specfit.parinfo)))
sp6e2abs.specfit(fittype='gaussian', guesses=[-1, 58, 1])
sp7e2abs.specfit(fittype='gaussian', guesses=[-1, 58, 1])
sp9e2abs.specfit(fittype='gaussian', guesses=[-1, 58, 1])
sp10e2abs.specfit(fittype='gaussian', guesses=[-1, 58, 1])
sp13e2abs.specfit(fittype='gaussian', guesses=[-1, 58, 1])
mlpi6e2abs = sp6e2abs.specfit.parinfo
mlpi7e2abs = sp7e2abs.specfit.parinfo
mlpi9e2abs = sp9e2abs.specfit.parinfo
mlpi10e2abs = sp9e2abs.specfit.parinfo
mlpi13e2abs = sp9e2abs.specfit.parinfo
print "gaussian 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp6e2abs.specfit.parinfo)))
print "gaussian 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp7e2abs.specfit.parinfo)))
print "gaussian 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp9e2abs.specfit.parinfo)))
print "gaussian 10-10 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp10e2abs.specfit.parinfo)))
print "gaussian 13-13 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp13e2abs.specfit.parinfo)))


"""
e2emi
hf 6-6 peak=0.034351759439 +/- 0.00118744160803, fwhm=5.68435483649 +/- 0.349980107582, integral=0.207855909422 +/- 0.0146764957049, center=56.4977738188 +/- 0.0667605561957
hf 7-7 peak=0.0166315025984 +/- 0.00130844782908, fwhm=5.62210537792 +/- 0.743037418582, integral=0.0995319882672 +/- 0.015308726158, center=56.4571801499 +/- 0.134505076546
hf 9-9 peak=0.00338848050734 +/- 0.000989794088329, fwhm=6.98406416945 +/- 1.93467593299, integral=0.0251910008648 +/- 0.0101411168447, center=56.4726962808 +/- 0.431893063304
gaussian 6-6 peak=0.0440648739987 +/- 0.000771763111177, fwhm=14.5271473445 +/- 0.293805143558, integral=0.68140463638 +/- 0.0182303648913, center=56.9990990941 +/- 0.124767557222
gaussian 7-7 peak=0.0457256919253 +/- 0.000881571520688, fwhm=10.8556466843 +/- 0.241669400766, integral=0.528382221133 +/- 0.0155608555916, center=56.297413005 +/- 0.102627547595
gaussian 9-9 peak=0.0228509869508 +/- 0.00077060378904, fwhm=11.5723391574 +/- 0.450629145718, integral=0.281486989149 +/- 0.0145002125858, center=56.2385609879 +/- 0.191364578003

e8
hf 6-6 peak=0.0552312979281 +/- 0.00339213398418, fwhm=12.0092504285 +/- 0.346442710066, integral=0.706046591242 +/- 0.0479084714681, center=59.8513485763 +/- 0.165718109489
hf 7-7 peak=0.0343580596483 +/- 0.00211133669018, fwhm=7.08292886544 +/- 0.477361456933, integral=0.259044108597 +/- 0.0236262598412, center=59.9800528242 +/- 0.291513634861
gaussian 6-6 peak=0.278613851186 +/- 0.00324112286397, fwhm=12.3888928639 +/- 0.166414956455, integral=3.67423906951 +/- 0.0652900918781, center=59.0587207674 +/- 0.0706699251669
gaussian 7-7 peak=0.19046626119 +/- 0.00375662817502, fwhm=11.1985631899 +/- 0.255043061619, integral=2.2704532914 +/- 0.0684040919618, center=59.7258654571 +/- 0.108306812353
gaussian 9-9 peak=0.0829717146021 +/- 0.00468072809892, fwhm=11.1981635835 +/- 0.729471923465, integral=0.989029128375 +/- 0.0852287264973, center=59.6499937385 +/- 0.309778206795

e2abs
hf 6-6 peak=-0.161810148064 +/- 0.000984062270758, fwhm=4.29488277478 +/- 0.0442752015202, integral=-0.739757335015 +/- 0.00885417492823, center=57.5484840785 +/- 0.0110282498426
hf 7-7 peak=-0.0642035748936 +/- 0.000934816176453, fwhm=4.51803495893 +/- 0.12446617443, integral=-0.308774201785 +/- 0.00962133785726, center=57.2170494137 +/- 0.0252379940539
hf 9-9 peak=-0.0097371317438 +/- 0.000271624308518, fwhm=6.49590535538 +/- 0.223255478522, integral=-0.0673291210331 +/- 0.00298031207001, center=57.3601759138 +/- 0.159023224602
gaussian 6-6 peak=-0.336181409247 +/- 0.00105372911297, fwhm=7.2878981115 +/- 0.0263769677346, integral=-2.60800365615 +/- 0.0124867865145, center=57.1430653294 +/- 0.0112012668288
gaussian 7-7 peak=-0.310241636013 +/- 0.000963459939972, fwhm=5.93209298899 +/- 0.0212721666744, integral=-1.95902619112 +/- 0.00929314446325, center=57.0861072678 +/- 0.00903345734588
gaussian 9-9 peak=-0.227156053636 +/- 0.00118632044585, fwhm=5.53519024599 +/- 0.0333797030406, integral=-1.33840970644 +/- 0.0106771857415, center=57.221370475 +/- 0.0141750545978
"""

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
flux_e2emi = {'66': (full_integral(hfpi6e2emi, mlpi6e2emi)*u.Jy).to(u.K, u.brightness_temperature(areas['e2emi'], freqs['66']))*u.km/u.s,
              '77': (full_integral(hfpi7e2emi, mlpi7e2emi)*u.Jy).to(u.K, u.brightness_temperature(areas['e2emi'], freqs['77']))*u.km/u.s,
              '99': (full_integral(hfpi9e2emi, mlpi9e2emi)*u.Jy).to(u.K, u.brightness_temperature(areas['e2emi'], freqs['99']))*u.km/u.s,
              '1010': (full_integral(None, mlpi10e2emi)*u.Jy).to(u.K, u.brightness_temperature(areas['e2emi'], freqs['1010']))*u.km/u.s,
              '1313': (full_integral(None, mlpi13e2emi)*u.Jy).to(u.K, u.brightness_temperature(areas['e2emi'], freqs['1313']))*u.km/u.s,
          }
taue2emi = {'66': 100,
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

Nu_e8 = {line: Nu_thin_hightex(flux_e8[line], linestrengths[line], freqs[line]) for line in freqs}
Nu_e2emi = {line: Nu_thin_hightex(flux_e2emi[line], linestrengths[line], freqs[line]) for line in freqs}
Nu_e8_taucorr = {line: Nu_thin_hightex(flux_e8[line], linestrengths[line], freqs[line], tau=taue8[line]) for line in freqs}
Nu_e2emi_taucorr = {line: Nu_thin_hightex(flux_e2emi[line], linestrengths[line], freqs[line], tau=taue2emi[line]) for line in freqs}

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

    print("Tex={0}, Ntot={1}, Q_rot={2}".format(tex, Ntot, Q_rot))

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
ax.plot(eups, np.exp(eups*slope + intcpt), color=a.get_color(), label='$T_{{ex}}={0:0.1f}$\n$N(\\mathrm{{NH}}_3)={1:0.1e}$'.format(tex, Ntot))
print("Mass lower limit: {0}".format((Ntot/XNH3 * areas['e8']*distance**2 * 2.8*u.Da).to(u.M_sun, u.dimensionless_angles())))
a,b,c = ax.errorbar(x=[eupper[line] for line in freqs],
            y=[Nu_e8_taucorr[line][0].value*degeneracy[line] for line in freqs],
            yerr=[Nu_e8_taucorr[line][1].value*degeneracy[line] for line in freqs],
            marker='s',
            linestyle='none')
Ntot,tex,slope,intcpt = fit_tex([eupper[line] for line in freqs],
                                [Nu_e8_taucorr[line][0].value*degeneracy[line] for line in freqs])
ax.plot(eups, np.exp(eups*slope + intcpt), color=a.get_color(), label='$T_{{ex}}={0:0.1f}$\n$N(\\mathrm{{NH}}_3)={1:0.1e}$'.format(tex, Ntot))
ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_ylabel("$N_u / g$ [cm$^{-2}$]")
ax.set_xlabel("$E_u$ [K]")
ax.set_title('e8')
pl.legend(loc='best')

print("e2emi")
fig = pl.figure(2)
fig.clf()
ax = pl.gca()
a,b,c = ax.errorbar(x=[eupper[line] for line in freqs],
            y=[Nu_e2emi[line][0].value*degeneracy[line] for line in freqs],
            yerr=[Nu_e2emi[line][1].value*degeneracy[line] for line in freqs],
            marker='s',
            linestyle='none')
Ntot,tex,slope,intcpt = fit_tex([eupper[line] for line in freqs], [Nu_e2emi[line][0].value*degeneracy[line] for line in freqs])
ax.plot(eups, np.exp(eups*slope + intcpt), color=a.get_color(), label='$T_{{ex}}={0:0.1f}$\n$N(\\mathrm{{NH}}_3)={1:0.1e}$'.format(tex, Ntot))
print("Mass lower limit: {0}".format((Ntot/XNH3 * areas['e2emi']*distance**2 * 2.8*u.Da).to(u.M_sun, u.dimensionless_angles())))
a,b,c = ax.errorbar(x=[eupper[line] for line in freqs],
            y=[Nu_e2emi_taucorr[line][0].value*degeneracy[line] for line in freqs],
            yerr=[Nu_e2emi_taucorr[line][1].value*degeneracy[line] for line in freqs],
            marker='s',
            linestyle='none')
Ntot,tex,slope,intcpt = fit_tex([eupper[line] for line in freqs],
                                [Nu_e2emi_taucorr[line][0].value*degeneracy[line] for line in freqs])
ax.plot(eups, np.exp(eups*slope + intcpt), color=a.get_color(), label='$T_{{ex}}={0:0.1f}$\n$N(\\mathrm{{NH}}_3)={1:0.1e}$'.format(tex, Ntot))
ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_ylabel("$N_u / g$ [cm$^{-2}$]")
ax.set_xlabel("$E_u$ [K]")
ax.set_title('e2emi')
pl.legend(loc='best')

pl.draw()
pl.show()


