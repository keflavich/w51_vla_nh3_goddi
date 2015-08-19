from astropy.io import fits
from astropy import wcs
from astropy import coordinates
import aplpy
import pylab as pl
import numpy as np

hdu = fits.open('W51-25GHzcont-pb.map.image.fits')
hdu[0].header = wcs.WCS(hdu[0].header).sub([wcs.WCSSUB_CELESTIAL]).to_header()
hdu[0].data = hdu[0].data.squeeze()
fig=pl.figure(1)
fig.clf()
F = aplpy.FITSFigure(hdu, figure=fig, subplot=[0.1,0.1, 0.7,0.8])
F.show_grayscale()
F.show_regions('../maserspots_circles.reg')

w51e2 = coordinates.SkyCoord.from_name('w51 e2')
F.recenter(w51e2.ra.deg, w51e2.dec.deg, radius=3/3600.)

F.show_contour('w51e2_340.CO_290chan.m0E.blue.mask.fits', colors=['b'],
               levels=np.linspace(5, 80, 5))
F.show_contour('w51e2_340.CO_290chan.m0E.red.mask.fits', colors=['r'],
               levels=np.linspace(5, 80, 5))

n,velocity,intensity,dra,ddec = np.loadtxt('../Etoka_masertable.txt').swapaxes(0,1)
norm = pl.mpl.colors.Normalize(vmin=velocity.min(), vmax=velocity.max())
ax = F._figure.add_axes([0.85, 0.1, 0.05, 0.8])
cb = pl.mpl.colorbar.ColorbarBase(ax, norm=norm, cmap=pl.cm.jet)
cb.set_label('$V_{LSR} ($km s$^{-1})$')
pl.draw()
pl.show()
