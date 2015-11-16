import numpy as np
import aplpy
from astropy.io import fits
from astropy import wcs

# Load the fitted parameter cube
# It has 14 slices: the first 7 are the fitted parameters, the next 7 are the
# errors on those parameters
# The parameters are:
# [ centroid velocity, hyperfine offset velocity 1, hyperfine amplitude 1,
#   hyperfine linewidth 1, hyperfine offset velocity 2, hyperfine amplitude 2,
#   hyperfine linewidth 2 ]
e2hf = fits.open('fitted_nh366_hf_emission.fits')

# Create a mask for the data: mask any points where the hyperfine amplitude is
# < 2.5 mJy (unless it is <-10 mJy, in which case it's absorption)
# and mask out any places where the line width is <0.1 km/s
e2hfd = e2hf[0].data
ok = (((e2hfd[2, :, :] > 0.0025) | (e2hfd[2, :, :] < -0.01))
      & (e2hfd[3, :, :] > 0.1)
      & ((e2hfd[5, :, :] > 0.0025) | (e2hfd[5, :, :] < -0.01))
      & (e2hfd[6, :, :] > 0.1))

e2hf[0].data[:, ~ok] = np.nan

# we have to fix the header; aplpy has problems using headers that have too
# many parameters in them
w = wcs.WCS(e2hf[0].header).sub([wcs.WCSSUB_CELESTIAL])
hdu = fits.PrimaryHDU(data=e2hf[0].data[0,:,:], header=w.to_header())
Fhf = aplpy.FITSFigure(hdu)
Fhf.show_colorscale(vmin=54, vmax=60)
Fhf.show_colorbar()
Fhf.save("e2_emi_velocity_field.eps")
Fhf.save("e2_emi_velocity_field.png")
Fhf.recenter(290.93295, 14.509594, radius=1/3600.)
Fhf.save("e2_abs_velocity_field.eps")
Fhf.save("e2_abs_velocity_field.png")

hdu = fits.PrimaryHDU(data=e2hf[0].data[3,:,:], header=w.to_header())
Fhfs1 = aplpy.FITSFigure(hdu)
Fhfs1.show_colorscale(vmin=0.4, vmax=4)
Fhfs1.show_colorbar()
Fhfs1.set_title('Width of inner HF lines')

hdu = fits.PrimaryHDU(data=e2hf[0].data[6,:,:], header=w.to_header())
Fhfs2 = aplpy.FITSFigure(hdu)
Fhfs2.show_colorscale(vmin=0.4, vmax=4)
Fhfs2.show_colorbar()
Fhfs2.set_title('Width of outer HF lines')




e2ml = fits.open('fitted_nh366_mainline_emission.fits')

e2mld = e2ml[0].data
okml = (((e2mld[0, :, :] > 0.0025) | (e2mld[0, :, :] < -0.01))
      & (e2mld[2, :, :] > 0.1))

e2ml[0].data[:, ~okml] = np.nan

w = wcs.WCS(e2ml[0].header).sub([wcs.WCSSUB_CELESTIAL])
hdu = fits.PrimaryHDU(data=e2ml[0].data[1,:,:], header=w.to_header())
Fml = aplpy.FITSFigure(hdu)
Fml.show_colorscale(vmin=54, vmax=60)
Fml.show_colorbar()





e8ml = fits.open('e8_fitted_nh366_mainline_emission.fits')

e8mld = e8ml[0].data
okml = (((e8mld[0, :, :] > 0.0025) | (e8mld[0, :, :] < -0.01))
      & (e8mld[2, :, :] > 0.1))

e8ml[0].data[:, ~okml] = np.nan

w = wcs.WCS(e8ml[0].header).sub([wcs.WCSSUB_CELESTIAL])
hdu = fits.PrimaryHDU(data=e8ml[0].data[1,:,:], header=w.to_header())
Fe8ml = aplpy.FITSFigure(hdu)
Fe8ml.show_colorscale(vmin=54, vmax=64)
Fe8ml.show_colorbar()



e8hf = fits.open('e8_fitted_nh366_hf_emission.fits')

e8hfd = e8hf[0].data
ok = (((e8hfd[2, :, :] > 0.0015) | (e8hfd[2, :, :] < -0.01))
      & (e8hfd[3, :, :] > 0.1)
      & ((e8hfd[5, :, :] > 0.0015) | (e8hfd[5, :, :] < -0.01))
      & (e8hfd[6, :, :] > 0.1))

e8hf[0].data[:, ~ok] = np.nan

w = wcs.WCS(e8hf[0].header).sub([wcs.WCSSUB_CELESTIAL])
hdu = fits.PrimaryHDU(data=e8hf[0].data[0,:,:], header=w.to_header())
Fe8hf = aplpy.FITSFigure(hdu)
Fe8hf.show_colorscale(vmin=54, vmax=64)
Fe8hf.show_colorbar()

hdu = fits.PrimaryHDU(data=e8hf[0].data[3,:,:], header=w.to_header())
Fe8hfs1 = aplpy.FITSFigure(hdu)
Fe8hfs1.show_colorscale(vmin=0.4, vmax=5)
Fe8hfs1.show_colorbar()
Fe8hfs1.set_title('Width of inner HF lines')

hdu = fits.PrimaryHDU(data=e8hf[0].data[6,:,:], header=w.to_header())
Fe8hfs2 = aplpy.FITSFigure(hdu)
Fe8hfs2.show_colorscale(vmin=0.4, vmax=5)
Fe8hfs2.show_colorbar()
Fe8hfs2.set_title('Width of outer HF lines')
