import os
from hf_only_model import hfonly_66_fixed_fitter, hfonly_fitter, sixsix_movinghf_fitter
from spectral_cube import SpectralCube
from astropy import units as u
import numpy as np
import pyspeckit

pyspeckit.fitters.default_Registry.add_fitter('hfonly', hfonly_fitter(), 7)
pyspeckit.fitters.default_Registry.add_fitter('hfonly66', hfonly_66_fixed_fitter(), 4)
pyspeckit.fitters.default_Registry.add_fitter('sixsix_movinghf', sixsix_movinghf_fitter(), 5)

cube = SpectralCube.read('/Volumes/passport/W51-GODDI/W51e2_66_baselined-sc-pb.cube.image.fits')
scube = cube[:, 230:307, 205:270]
errmap = scube.spectral_slab(-20*u.km/u.s, 10*u.km/u.s).std(axis=0)
mn = scube.min(axis=0).value


guesses = np.empty((7, scube.shape[1], scube.shape[2]))
guesses[0,:,:] = 58
guesses[1,:,:] = 26.9
guesses[2,:,:] = 0.010
guesses[3,:,:] = 1.5
guesses[4,:,:] = 31.4
guesses[5,:,:] = 0.010
guesses[6,:,:] = 1.5

negmask = mn<-0.005
guesses[2, negmask] = -0.1
guesses[5, negmask] = -0.1


#peak = scube.max(axis=0)
#mask = peak>0.005*peak.unit

pcube = pyspeckit.Cube(cube=scube)
parfilename = 'fitted_nh366_hf_emission.fits'
if os.path.exists(parfilename):
    pcube.load_model_fit(parfilename, 7, 1, 'hfonly')
else:
    pcube.fiteach(fittype='hfonly',
                  guesses=guesses,
                  limits=[(50, 70), (24, 29), (-10, 5), (0, 4), (29, 34), (-10, 5), (0, 4)],
                  limitedmin=[True]*7,
                  limitedmax=[True]*7,
                  integral=False, errmap=errmap.value,
                  signal_cut=0,
                  #maskmap=mask,
                  multicore=3)
    pcube.write_fit(parfilename, clobber=True)
    #pcube.fiteach(fittype='hfonly66', guesses=[58, 1, 1, 1],
    #              integral=False, errmap=errmap.value, maskmap=mask, multicore=3)

gguesses = np.empty((3, scube.shape[1], scube.shape[2]))
gguesses[0,:,:] = 0.010
gguesses[0,negmask] = -0.10
gguesses[1,:,:] = 58
gguesses[2,:,:] = 1.5

pcube2 = pyspeckit.Cube(cube=scube)
parfilename = 'fitted_nh366_mainline_emission.fits'
if os.path.exists(parfilename):
    pcube2.load_model_fit(parfilename, 3, 1, 'gaussian')
else:
    pcube2.fiteach(fittype='gaussian', guesses=gguesses,
                   limits=[(-10, 5), (50, 70), (0,5)],
                   limitedmin=[True]*3,
                   limitedmax=[True]*3,
                   integral=False,
                   errmap=errmap.value,
                   signal_cut=0,
                   #maskmap=mask,
                   multicore=3)
    pcube2.write_fit(parfilename, clobber=True)

scubee8 = cube[:, 40:140, 200:300]
errmape8 = scubee8.spectral_slab(-20*u.km/u.s, 10*u.km/u.s).std(axis=0)

pcube_e8_2 = pyspeckit.Cube(cube=scubee8)
parfilename = 'e8_fitted_nh366_mainline_emission.fits'
if os.path.exists(parfilename):
    pcube_e8_2.load_model_fit(parfilename, 3, 1, 'gaussian')
else:
    pcube_e8_2.fiteach(fittype='gaussian', guesses=[1, 58, 1],
                   limits=[(-10, 5), (50, 70), (0,7)],
                   limitedmin=[True]*3,
                   limitedmax=[True]*3,
                   integral=False,
                   errmap=errmape8.value,
                   signal_cut=0,
                   #maskmap=mask,
                   multicore=3)
    pcube_e8_2.write_fit(parfilename, clobber=True)

pcube_e8 = pyspeckit.Cube(cube=scubee8)
parfilename = 'e8_fitted_nh366_hf_emission.fits'
if os.path.exists(parfilename):
    pcube_e8_2.load_model_fit(parfilename, 7, 1, 'hfonly')
else:
    pcube_e8.fiteach(fittype='hfonly', guesses=[58, 26, 1, 2, 32, 1, 2],
                  limits=[(50, 70), (24, 29), (-10, 5), (0, 7), (29, 34), (-10, 5), (0, 7)],
                  limitedmin=[True]*7,
                  limitedmax=[True]*7,
                  integral=False, errmap=errmape8.value,
                   signal_cut=0,
                  #maskmap=mask,
                  multicore=3)
    pcube_e8.write_fit(parfilename, clobber=True)
    #pcube_e8.fiteach(fittype='hfonly66', guesses=[58, 1, 1, 1],
    #              integral=False, errmap=errmap.value, maskmap=mask, multicore=3)
