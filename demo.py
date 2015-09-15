import numpy, sys
import matplotlib.pyplot as plt
import MINERVAphot as mp

# ---------------------#
# First calibrate a defined list of files using lists of calibration files
# ---------------------#
biaslist = numpy.loadtxt('biases.list', dtype='object')
darklist = numpy.loadtxt('darks.list', dtype='object')
flatlist = numpy.loadtxt('flats.list', dtype='object')
unproclist = numpy.loadtxt('images.list', dtype='object')
Saveloc = "./CalibList/"

#procimgs = mp.CalibrateList(unproclist, biaslist, darklist, flatlist, Saveloc)
#print procimgs


# ---------------------#
# Calibrate a night of observations based on a target file.
# This specific example presumes that there is a directory named "n20150411.T3"
# and a target file named "n20150411.T3.txt" present in the location the program is run.
# Right now it assumes that the directory and the target filename are the same.
# ---------------------#
night = "n20150411.T3"
telno = "T3"
procimgs, filters = mp.CalibrateNight(night, telno, skipbias=True, skipdark=True, skipflat=True, skipproc=True)
# ---------------------#
# This returns a dictionary with lists of the processed image locations
# and a list of all the filters it detected
# ---------------------#

print procimgs['ip']
print filters
print procimgs[filters[0]]

KELT3gp = procimgs[filters[0]] # specific to the night of n20150411

AperProperties = {}
AperProperties['L'] = 10 # Half-size of the box for re-centering, in pixels
AperProperties['rAp'] = 12 # radius of the phot aperture, in pixels
AperProperties['rBkgIn'] = 15 # inner radius of the bkg annulus, in pixels
AperProperties['rBkgOut'] = 20 # outer radius of the bkg annulus, in pixes

StarLocs = mp.GraphicalPickLocs(KELT3gp[0])
AperProperties['rastar'] = StarLocs[0]/15. # to decimal degress
AperProperties['decstar'] = StarLocs[1]

if AperProperties['rastar'].size > 1:
	photresults = mp.MultiAper(KELT3gp, AperProperties)
else:
	photresults = mp.SingleAper(KELT3gp, AperProperties)

time = photresults[:,0]
detrended, error = mp.SimpleDetrend(photresults) # this assumes the first aperture is the target, and uses all the others as comparison stars

plt.errorbar(time, detrended, yerr=error, fmt='.k', capsize=0)
plt.show()

