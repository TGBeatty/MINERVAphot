import numpy, sys, aplpy, astropy
import matplotlib.pyplot as plt
import MINERVAphot as mp

def on_click(event):
    if event.inaxes is not None:
    	WCSloc = w.wcs_pix2world(event.xdata, event.ydata, 0)
    	fig.show_markers(WCSloc[0], WCSloc[1])
    	fig.refresh()
    else:
        print 'Clicked ouside axes bounds but inside plot window'

AperProperties = {}
AperProperties['rastar'] = 23.2329888 # decimal hours
AperProperties['decstar'] = 8.76127777 # decimal deg.
AperProperties['L'] = 5 # Half-size of the box for centering
AperProperties['rAp'] = 12 # radius of the phot aperture
AperProperties['rBkgIn'] = 15 # inner radius of the bkg annulus
AperProperties['rBkgOut'] = 20 # outer radius of the bkg annulus

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

photresults = mp.SingleAper(proclist, AperProperties)

#photresults = mp.MultiAper(proclist, AperProperties)

print photresults

numpy.savetxt('wasp52.lc', photresults)

plt.plot(photresults[:,0], photresults[:,1], '.')
plt.show()