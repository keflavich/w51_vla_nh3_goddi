from pyspeckit.spectrum.models.inherited_gaussfitter import gaussian
from pyspeckit.spectrum.models import model
from astropy import units as u
import numpy as np

def hfonly(x, dv, dx1, amp1, width1, dx2, amp2, width2, **kwargs):

    return (gaussian(x, amp1, dv+dx1, width1)+
            gaussian(x, amp1, dv-dx1, width1)+
            gaussian(x, amp2, dv+dx2, width2)+
            gaussian(x, amp2, dv-dx2, width2))

def hfonly_fitter():

    myclass = model.SpectralModel(hfonly, 7,
            parnames=['shift',
                      'hf1offset', 'hf1amp', 'hf1width',
                      'hf2offset', 'hf2amp', 'hf2width',
                     ],
            parlimited=[(False, False),
                        (True, False),(False,False),(True,False),
                        (True, False),(False,False),(True,False),
                       ],
            parlimits=[(0,0), (0,0), (0,0), (0,0), (0,0), (0,0), (0,0), ],
            shortvarnames=(r'\Delta x',
                           r'\Delta hf_1', r'A_1', r'\sigma_1',
                           r'\Delta hf_2', r'A_2', r'\sigma_2',
                          ),
            centroid_par='shift',
           )
    myclass.__name__ = "hfonly"

    return myclass

def hfonly_66_fixed(x, dv, amp1, amp2, width, **kwargs):

    dx1 = 26.9
    dx2 = 31.4

    return (gaussian(x, amp1, dv+dx1, width)+
            gaussian(x, amp1, dv-dx1, width)+
            gaussian(x, amp2, dv+dx2, width)+
            gaussian(x, amp2, dv-dx2, width))

def hfonly_66_fixed_fitter():

    myclass = model.SpectralModel(hfonly_66_fixed, 4,
            parnames=['shift', 'hf1amp', 'hf2amp', 'width', ],
            parlimited=[(False, False), (False,False), (False,False),
                        (True,False),
                       ],
            parlimits=[(0,0), (0,0), (0,0), (0,0), ],
            shortvarnames=(r'\Delta x', r'A_1', r'A_2', r'\sigma',
                          ),
            centroid_par='shift',
           )
    myclass.__name__ = "hfonly_66_fixed"

    return myclass

def sixsix_movinghf(xarr, dv, tau, dx1, dx2, width, tex=400,
                    tbackground=8200, strength1=0.0081, strength2=0.0081,
                    **kwargs):
    """
    Try to fit the optical depth by allowing the hyperfine lines to move a little
    """

    xarr = xarr.to(u.km/u.s).value
    tau_nu_cumul = np.zeros_like(xarr)

    for offset,width,strength in ((dx1, width, strength1),
                                  (dx2, width, strength2),
                                  (0, width, 1)):


        if offset > 0:
            frac = (np.array(strength *
                               np.exp(-(xarr-dv+offset)**2 /
                                       (2.0*width**2)))
                      + np.array(strength *
                                 np.exp(-(xarr-dv-offset)**2 /
                                        (2.0*width**2)))
                     )
        else:
            frac = np.array(strength *
                              np.exp(-(xarr-dv+offset)**2 /
                                      (2.0*width**2)))
        frac[frac!=frac] = 0 # avoid nans

        tau_nu_cumul += frac

    spec = (1.0-np.exp(-np.array(tau_nu_cumul*tau)))*(tex-tbackground)

    return spec

def sixsix_movinghf_fitter():

    myclass = model.SpectralModel(sixsix_movinghf, 7,
            parnames=['shift', 'tau', 'hf1off', 'hf2off', 'width',
                      'tex', 'tbackground', ],
            parlimited=[(False, False), (True,False), (True,False),
                        (True,False), (True,False),
                        (True,False), (True, False),
                       ],
            parlimits=[(0,0), (0,0), (0,0), (0,0), (0,0), (0,0), (0,0)],
            shortvarnames=(r'\Delta x', r'\tau', r'\Delta_{hf1}',
                           r'\Delta_{hf2}', r'\sigma',
                           r'T_{ex}', r'T_{bg}',
                          ),
            centroid_par='shift',
           )
    myclass.__name__ = "sixsix_movinghf"

    return myclass

def sevenseven_movinghf(xarr, dv, tau, dx1, dx2, width, tex=400,
                    tbackground=8200, strength1=0.0060, strength2=0.0060,
                    **kwargs):
    """
    Try to fit the optical depth by allowing the hyperfine lines to move a little
    """

    xarr = xarr.to(u.km/u.s).value
    tau_nu_cumul = np.zeros_like(xarr)

    for offset,width,strength in ((dx1, width, strength1),
                                  (dx2, width, strength2),
                                  (0, width, 1)):


        if offset > 0:
            frac = (np.array(strength *
                               np.exp(-(xarr-dv+offset)**2 /
                                       (2.0*width**2)))
                      + np.array(strength *
                                 np.exp(-(xarr-dv-offset)**2 /
                                        (2.0*width**2)))
                     )
        else:
            frac = np.array(strength *
                              np.exp(-(xarr-dv+offset)**2 /
                                      (2.0*width**2)))
        frac[frac!=frac] = 0 # avoid nans

        tau_nu_cumul += frac

    spec = (1.0-np.exp(-np.array(tau_nu_cumul*tau)))*(tex-tbackground)

    return spec

def sevenseven_movinghf_fitter():

    myclass = model.SpectralModel(sevenseven_movinghf, 7,
            parnames=['shift', 'tau', 'hf1off', 'hf2off', 'width',
                      'tex', 'tbackground', ],
            parlimited=[(False, False), (True,False), (True,False),
                        (True,False), (True,False),
                        (True,False), (True, False),
                       ],
            parlimits=[(0,0), (0,0), (0,0), (0,0), (0,0), (0,0), (0,0)],
            shortvarnames=(r'\Delta x', r'\tau', r'\Delta_{hf1}',
                           r'\Delta_{hf2}', r'\sigma',
                           r'T_{ex}', r'T_{bg}',
                          ),
            centroid_par='shift',
           )
    myclass.__name__ = "sevenseven_movinghf"

    return myclass


def ninenine_hf(xarr, dv, tex, tau, width, dx1=27.0, dx2=30.1,
                tbackground=0, strength1=0.0037, strength2=0.0037,
                **kwargs):

    xarr = xarr.to(u.km/u.s).value
    tau_nu_cumul = np.zeros_like(xarr)

    for offset,width,strength in ((dx1, width, strength1),
                                  (dx2, width, strength2),
                                  (0, width, 1)):


        if offset > 0:
            frac = (np.array(strength *
                               np.exp(-(xarr-dv+offset)**2 /
                                       (2.0*width**2)))
                      + np.array(strength *
                                 np.exp(-(xarr-dv-offset)**2 /
                                        (2.0*width**2)))
                     )
        else:
            frac = np.array(strength *
                              np.exp(-(xarr-dv+offset)**2 /
                                      (2.0*width**2)))
        frac[frac!=frac] = 0 # avoid nans

        tau_nu_cumul += frac

    spec = (1.0-np.exp(-np.array(tau_nu_cumul*tau)))*(tex-tbackground)

    return spec

def ninenine_hf_fitter():

    myclass = model.SpectralModel(ninenine_hf, 4,
            parnames=['shift', 'tex', 'tau', 'width' ],
            parlimited=[(False, False), (True,False), (True,False),
                        (True,False)
                       ],
            parlimits=[(0,0), (0,0), (0,0), (0,0),],
            shortvarnames=(r'\Delta x', r'T_{ex}', r'\tau',
                           r'\sigma',
                          ),
            centroid_par='shift',
           )
    myclass.__name__ = "ninenine_hf"

    return myclass
