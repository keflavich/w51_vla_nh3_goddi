## Code README ##

``fit_hyperfines_and_report.py``

Fit the hyperfine lines and gaussian main component to each spectrum.
Hyperfines are fit independently from the main line.  Also does rotation
(Boltzmann) diagram fitting and mass estimation for the *emission* sources.
Relies on `goddi_nh3_fits.py`, [pyspeckit](pyspeckit.bitbucket.org), `hf_only_model.py`,
and [astroquery](pyspeckit.bitbucket.org).

``goddi_nh3_fits.py``

Loads the spectra and converts them to Kelvins.  Can do some sophisticated NH3 fitting.
Most of the work in the `__main__` function isn't used any more, though - we just use
this code for reading the data.

``etoka_overlay.py``

Overplot the CH3OH masers


``pubfiguresrc``

fancy plotting setup (thicker lines, bigger text)


``nh3_abs_trot.py``

The absorption-based rotation diagram fitting.  Has hard-coded fluxes!


``hf_only_model.py``

A variatety of pyspeckit models for fitting the NH3 hyperfines

``fit_66_cube.py``
``plot_66_maps.py``

Fit the 6-6 lines pixel-by-pixel and then create velocity field maps
