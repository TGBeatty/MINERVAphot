import numpy
import matplotlib.pyplot as plt
import MINERVAphot as mp

rphotap = 12
rbkgin = 15
rbkgout = 20

biaslist = numpy.loadtxt('biases.list', dtype='object')
darklist = numpy.loadtxt('darks.list', dtype='object')
flatlist = numpy.loadtxt('flats.list', dtype='object')
unproclist = numpy.loadtxt('unreduced.list', dtype='object')

biasloc = mp.MakeBias(biaslist, './CalibFiles/')
darkloc = mp.MakeDark(darklist, './CalibFiles/', biasloc)
flatloc = mp.MakeFlat(flatlist, './CalibFiles/', biasloc, darkloc)
proclist = mp.BiasDarkFlat(unproclist, './Processed/', biasloc, darkloc, flatloc)

photresults = mp.SingleAper(proclist, rphotap, rbkgin, rbkgout)

numpy.savetxt('wasp52.lc', photresults)

#fitslist = numpy.loadtxt('wasp52.list', dtype='object')
#OGphotresults = minphot.SingleAper(fitslist, rphotap, rbkgin, rbkgout)

plt.plot(photresults[:,0], photresults[:,1], '.')
plt.show()