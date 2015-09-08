# MINERVAphot

This is a python-based photometry pipeline for MINERVA.

Right now, it can do:

---------------------------------------

1) Create a master bias image from a list of input files.

2) Create a master dark image from a list of input files.

3) Create a master flat image from a list of input files.

4) Bias/Dark/Flat subtract a list of science images.

5) Allows the user to graphically select the location of the photometric apertures, and centers automatically.

6) Perform single aperture photometry if only one aperture is specified.

6) Perform multi aperture photometry if more than one aperture is specified.

---------------------------------------

Otherwise, it requires the user to specify the location of the target star(s), the size of the photometric aperture, and the size of the background annulus. Right now it uses the same aperture and annulus sizes for everything.

The output is a numpy array of the form:

time, flux, error, airmass, x-center, y-center

For more than one star, the 'flux', 'error', and 'center' columns will expand to match the number of number of stars.
