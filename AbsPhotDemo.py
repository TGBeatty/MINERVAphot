import math, numpy, sys, pandas
import matplotlib.pyplot as plt
import MINERVAphot as mp
import statsmodels.formula.api as sm

targgplist = numpy.loadtxt('target.gp.list', unpack=True, dtype='object')
targrplist = numpy.loadtxt('target.rp.list', unpack=True, dtype='object')
targiplist = numpy.loadtxt('target.ip.list', unpack=True, dtype='object')
gplist = numpy.loadtxt('GD71.gp.list', unpack=True, dtype='object')
rplist = numpy.loadtxt('GD71.rp.list', unpack=True, dtype='object')
iplist = numpy.loadtxt('GD71.ip.list', unpack=True, dtype='object')
name, ra, dec, up, sigup, nup, gp, siggp, ngp, rp, sigrp, nrp, ip, sigip, nip, zp, sigzp, nzp = numpy.loadtxt('GD71.txt.clean.v3', unpack=True, skiprows=1)

AperProperties = {}
AperProperties['L'] = 5 # Half-size of the box for re-centering, in pixels
AperProperties['rAp'] = 15 # radius of the phot aperture, in pixels
AperProperties['rBkgIn'] = 20 # inner radius of the bkg annulus, in pixels
AperProperties['rBkgOut'] = 30 # outer radius of the bkg annulus, in pixes

AperProperties['rastar'] = ra[0:10]/15. # to decimal hours
AperProperties['decstar'] = dec[0:10]

gpphot = mp.MultiAper(gplist, AperProperties)
rpphot = mp.MultiAper(rplist, AperProperties)
ipphot = mp.MultiAper(iplist, AperProperties)

#numpy.savetxt('photresults.gp.dat', gpphot)
#numpy.savetxt('photresults.rp.dat', rpphot)
#numpy.savetxt('photresults.ip.dat', ipphot)
#gpphot = numpy.loadtxt('photresults.gp.dat')
#rpphot = numpy.loadtxt('photresults.rp.dat')
#ipphot = numpy.loadtxt('photresults.ip.dat')

gpmag = gp[0:10]
gpmagerr = siggp[0:10]
rpmag = rp[0:10]
rpmagerr = sigrp[0:10]
ipmag = ip[0:10]
ipmagerr = sigip[0:10]

AperProperties['rastar'] = 6.992108 # to decimal hours
AperProperties['decstar'] = -4.0910666

targgp = mp.SingleAper(targgplist, AperProperties)
targrp = mp.SingleAper(targrplist, AperProperties)
targip = mp.SingleAper(targiplist, AperProperties)

#numpy.savetxt('target.gp.dat', targgp)
#numpy.savetxt('target.rp.dat', targrp)
#numpy.savetxt('target.ip.dat', targip)
#targgp = numpy.loadtxt('target.gp.dat')
#targrp = numpy.loadtxt('target.rp.dat')
#targip = numpy.loadtxt('target.ip.dat')

target_gp, target_gp_err = mp.SingleAbsPhot(targgp, gpphot, gpmag, gpmagerr, targrp, rpphot, rpmag, rpmagerr)
print target_gp, target_gp_err

target_rp, target_rp_err = mp.SingleAbsPhot(targrp, rpphot, rpmag, rpmagerr, targgp, gpphot, gpmag, gpmagerr)
print target_rp, target_rp_err

target_ip, target_ip_err = mp.SingleAbsPhot(targip, ipphot, ipmag, ipmagerr, targgp, gpphot, gpmag, gpmagerr)
print target_ip, target_ip_err





