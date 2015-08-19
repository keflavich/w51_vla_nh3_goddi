"""
Code for populating table 3
"""
import numpy as np
import pyspeckit
import goddi_nh3_fits
from hf_only_model import hfonly_fitter


sp6e2emi,spK6e2emi,_,_ = goddi_nh3_fits.load_spectrum(6, object='w51e2-proto', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp6e2emi.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp7e2emi,spK7e2emi,_,_ = goddi_nh3_fits.load_spectrum(7, object='w51e2-proto', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp7e2emi.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp9e2emi,spK9e2emi,_,_ = goddi_nh3_fits.load_spectrum(9, object='w51e2-proto', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp9e2emi.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)

F = False
T = True
sp6e2emi.specfit(fittype='hfonly', guesses=[58, 26.9, 0.02, 2.0, 31.4, 0.02, 2.0], fixed=[F,F,F,F,F,F,F])
sp7e2emi.specfit(fittype='hfonly', guesses=[58, 26.9, 0.02, 2.0, 31.4, 0.02, 2.0], fixed=[F,F,F,F,F,F,F])
sp9e2emi.specfit(fittype='hfonly', guesses=[58, 27.0, 0.02, 2.0, 30.1, 0.02, 2.0], fixed=[F,T,F,F,T,F,F])

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

print "e2emi"
print "hf 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp6e2emi.specfit.parinfo)))
print "hf 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp7e2emi.specfit.parinfo)))
print "hf 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp9e2emi.specfit.parinfo)))
sp6e2emi.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp7e2emi.specfit(fittype='gaussian', guesses=[1, 58, 1])
sp9e2emi.specfit(fittype='gaussian', guesses=[1, 58, 1])
print "gaussian 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp6e2emi.specfit.parinfo)))
print "gaussian 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp7e2emi.specfit.parinfo)))
print "gaussian 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp9e2emi.specfit.parinfo)))

sp6e8,spK6e8,_,_ = goddi_nh3_fits.load_spectrum(6, object='w51e8', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp6e8.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp7e8,spK7e8,_,_ = goddi_nh3_fits.load_spectrum(7, object='w51e8', headerfile='/Users/adam/work/w51/goddi/W51-25GHzcont.map.image.fits')
sp7e8.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)
sp9e8,spK9e8,_,_ = goddi_nh3_fits.load_spectrum(9, object='w51e8', headerfile='/Users/adam/work/w51/goddi/W51-27GHzcont.map.image.fits')
sp9e8.specfit.Registry.add_fitter('hfonly',hfonly_fitter(),7)

sp6e8.plotter()
sp7e8.plotter()
sp9e8.plotter()
limits = [(50, 65), (25,30), (0, 1), (0.1, 6), (30, 35), (0, 1), (0.1, 6)]
limited = [(True,True)]*7
tied = ['']*7
tied[6] = 'p[3]'
tied[5] = 'p[2]'
sp6e8.specfit(fittype='hfonly', guesses=[58, 26.9, 0.02, 2.0, 31.4, 0.02, 2.0], fixed=[F,F,F,F,F,F,F], limited=limited, limits=limits)
sp7e8.specfit(fittype='hfonly', guesses=[58, 26.9, 0.02, 2.0, 31.4, 0.02, 2.0], fixed=[F,F,F,F,F,F,F], limited=limited, limits=limits, tied=tied)
sp9e8.specfit(fittype='hfonly', guesses=[58, 27.0, 0.02, 2.0, 30.1, 0.02, 2.0], fixed=[F,T,F,F,T,F,F], limited=limited, limits=limits, tied=tied)
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
print "gaussian 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp6e8.specfit.parinfo)))
print "gaussian 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp7e8.specfit.parinfo)))
print "gaussian 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp9e8.specfit.parinfo)))
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
print "e2abs"
print "hf 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp6e2abs.specfit.parinfo)))
print "hf 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp7e2abs.specfit.parinfo)))
print "hf 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals(sp9e2abs.specfit.parinfo)))
sp6e2abs.specfit(fittype='gaussian', guesses=[-1, 58, 1])
sp7e2abs.specfit(fittype='gaussian', guesses=[-1, 58, 1])
sp9e2abs.specfit(fittype='gaussian', guesses=[-1, 58, 1])
print "gaussian 6-6 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp6e2abs.specfit.parinfo)))
print "gaussian 7-7 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp7e2abs.specfit.parinfo)))
print "gaussian 9-9 peak={0} +/- {1}, fwhm={2} +/- {3}, integral={4} +/- {5}, center={6} +/- {7}".format(*(tbl_vals_gaussian(sp9e2abs.specfit.parinfo)))


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
