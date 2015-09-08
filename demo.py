import numpy, sys, aplpy, astropy
import matplotlib.pyplot as plt
import MINERVAphot as mp

AperProperties = {}
AperProperties['rastar'] = 23.2329888 # decimal hours
AperProperties['decstar'] = 8.76127777 # decimal deg.
AperProperties['L'] = 10 # Half-size of the box for re-centering, in pixels
AperProperties['rAp'] = 12 # radius of the phot aperture, in pixels
AperProperties['rBkgIn'] = 15 # inner radius of the bkg annulus, in pixels
AperProperties['rBkgOut'] = 20 # outer radius of the bkg annulus, in pixes

biaslist = numpy.loadtxt('biases.list', dtype='object')
darklist = numpy.loadtxt('darks.list', dtype='object')
flatlist = numpy.loadtxt('flats.list', dtype='object')
unproclist = numpy.loadtxt('unreduced.list', dtype='object')

CalibFileSaveLoc = './CalibFiles/'

biasloc = mp.MakeBias(biaslist, CalibFileSaveLoc, skip=True)
darkloc = mp.MakeDark(darklist, CalibFileSaveLoc, biasloc, skip=True)
flatloc = mp.MakeFlat(flatlist, CalibFileSaveLoc, biasloc, darkloc, skip=True)
proclist = mp.BiasDarkFlat(unproclist, './Processed/', biasloc, darkloc, flatloc, skip=True)

StarLocs = mp.GraphicalPickLocs(proclist[0])
AperProperties['rastar'] = StarLocs[0]/15. # to decimal degress
AperProperties['decstar'] = StarLocs[1]

if AperProperties['rastar'].size > 1:
	photresults = mp.MultiAper(proclist, AperProperties)
else:
	photresults = mp.SingleAper(proclist, AperProperties)

numpy.savetxt('wasp52.lc', photresults)

plt.plot(photresults[:,0], photresults[:,1]/photresults[:,2], '.b')
plt.show()