# MINERVAphot

Last update: 09/22/15, v0.5

This is a python-based photometry pipeline for MINERVA.

Right now, it can do:

---------------------------------------

1) Calibrate a list of images based on a list of calibration (bias/dark/flat) images.

2) Calibrate a night's worth of images based on a MINERVA target file.

3) Do single and multiple aperture photometry.

4) Allows the user to graphically select the location of the photometric apertures, and centers automatically.

5) Does a simple detrending of mult-aperture results.

6) Do absolute photometry of a single target star using a set of Sloan standard stars. See AbsPhotDemo.py for an example.

---------------------------------------

Otherwise, it requires the user to specify the location of the target star(s), the size of the photometric aperture, and the size of the background annulus. Right now it uses the same aperture and annulus sizes for everything.

The aperture photometry output is a numpy array of the form:

time, airmass, flux, error, x-center, y-center

For more than one star, the 'flux', 'error', and 'center' columns will expand to match the number of number of stars. That is, 3 apertures will give you 3 columns of flux, followed by 3 columns of errors, etc.
