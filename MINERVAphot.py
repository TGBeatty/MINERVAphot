import numpy, progressbar, astropy, math, sys, os
import photutils as pu

def MakeBias(fitslist, outloc):
	if not os.path.exists(outloc):
		os.makedirs(outloc)

	nimages = len(fitslist)

	print "Making bias image..."
	bar = progressbar.ProgressBar(maxval=nimages, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	for i in range(nimages):
		hdulist = astropy.io.fits.open(fitslist[i])
		if i == 0:
			size = (nimages,hdulist[0].data.shape[0],hdulist[0].data.shape[1])
			tomedian = numpy.zeros(size)
		tomedian[i,:,:] = hdulist[0].data
		bar.update(i+1)
	bar.finish()

	hdulist[0].data = numpy.median(tomedian, axis=0)
	biasloc = outloc+"mbias.fits"
	hdulist.writeto(biasloc, clobber=True)

	return biasloc

def MakeDark(fitslist, outloc, mbiasloc):
	if not os.path.exists(outloc):
		os.makedirs(outloc)

	nimages = len(fitslist)
	hdulist = astropy.io.fits.open(mbiasloc)
	bias = hdulist[0].data

	print "Making dark image..."
	bar = progressbar.ProgressBar(maxval=nimages, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	for i in range(nimages):
		hdulist = astropy.io.fits.open(fitslist[i])
		if i == 0:
			size = (nimages,hdulist[0].data.shape[0],hdulist[0].data.shape[1])
			tomedian = numpy.zeros(size)
		tomedian[i,:,:] = hdulist[0].data - bias
		bar.update(i+1)
	bar.finish()

	hdulist[0].data = numpy.median(tomedian, axis=0)
	darkloc = outloc+"mdark.fits"
	hdulist.writeto(darkloc, clobber=True)

	return darkloc

def MakeFlat(fitslist, outloc, mbiasloc, mdarkloc):
	if not os.path.exists(outloc):
		os.makedirs(outloc)

	nimages = len(fitslist)
	hdulist = astropy.io.fits.open(mbiasloc)
	bias = hdulist[0].data
	hdulist = astropy.io.fits.open(mdarkloc)
	dark = hdulist[0].data

	print "Making flat image..."
	bar = progressbar.ProgressBar(maxval=nimages, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	for i in range(nimages):
		hdulist = astropy.io.fits.open(fitslist[i])
		if i == 0:
			size = (nimages,hdulist[0].data.shape[0],hdulist[0].data.shape[1])
			tomedian = numpy.zeros(size)
		tomedian[i,:,:] = (hdulist[0].data - bias - dark) / numpy.median(hdulist[0].data)
		bar.update(i+1)
	bar.finish()

	hdulist[0].data = numpy.median(tomedian, axis=0)
	flatloc = outloc+"mflat.fits"
	hdulist.writeto(flatloc, clobber=True)

	return flatloc

def BiasDarkFlat(fitslist, outloc, mbiasloc, mdarkloc, mflatloc):
	if not os.path.exists(outloc):
		os.makedirs(outloc)

	nimages = len(fitslist)
	procimgloc = numpy.empty(nimages, dtype='object') # What will be the list of processed images

	hdulist = astropy.io.fits.open(mbiasloc)
	bias = hdulist[0].data
	hdulist = astropy.io.fits.open(mdarkloc)
	dark = hdulist[0].data
	hdulist = astropy.io.fits.open(mflatloc)
	flat = hdulist[0].data 

	print "Calibrating images..."
	bar = progressbar.ProgressBar(maxval=nimages, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	for i in range(nimages):
		hdulist = astropy.io.fits.open(fitslist[i])
		image = hdulist[0].data
		hdulist[0].data = (image - bias - dark) / flat
		splitname = fitslist[i].split('/')
		filenamesplit = splitname[-1].partition('.fit')
		procimgloc[i] = outloc+filenamesplit[0]+".proc"+filenamesplit[1]+filenamesplit[2]
		hdulist.writeto(procimgloc[i], clobber=True)
		bar.update(i+1)
	bar.finish()

	return procimgloc

def HowellCenter(image, center, L):
	# Howell centroiding
	# L is the half-size of the subframe extracted to do the centering
	subimage = image[center[1]-2*L:center[1]+2*L,center[0]-2*L:center[0]+2*L]

	xpixels = numpy.arange(subimage.shape[1])
	ypixels = numpy.arange(subimage.shape[0])
	I = numpy.sum(subimage, axis=0)
	J = numpy.sum(subimage, axis=1)

	Isub = I-numpy.sum(I)/I.size
	Isub[Isub<0] = 0
	Jsub = J-numpy.sum(J)/J.size
	Jsub[Jsub<0] = 0
	xc = numpy.sum(Isub*xpixels)/numpy.sum(Isub)+center[0]-2*L
	yc = numpy.sum(Jsub*ypixels)/numpy.sum(Jsub)+center[1]-2*L
	centerloc = (xc, yc)

	return centerloc

def SingleAper(fitslist, rphot, rbkgin, rbkgout):
	rastar = 348.49483333
	decstar = 8.76127777
	L = 5. # centering box half size

	nimages = len(fitslist)
	jdutc = numpy.zeros(nimages)
	flux = numpy.zeros(nimages)
	error = numpy.zeros(nimages)
	xcenter = numpy.zeros(nimages)
	ycenter = numpy.zeros(nimages)
	airmass = numpy.zeros(nimages)

	print "Extracting photometry..."
	bar = progressbar.ProgressBar(maxval=nimages, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	for i in range(nimages):
		hdulist = astropy.io.fits.open(fitslist[i])
		image = hdulist[0].data 
		starttime = numpy.float(hdulist[0].header['JD'])
		exptime = numpy.float(hdulist[0].header['EXPTIME'])
		jdutc[i] = starttime + ((exptime/2.)/(24.*60*60.))
		alt = numpy.float(hdulist[0].header['OBJCTALT'])*math.pi/180.
		airmass[i] = 1./numpy.cos((math.pi/2.)-alt)

		try:
			wcscheck = hdulist[0].header['WCSAXES']
			w = astropy.wcs.WCS(hdulist[0].header)
			WCScenter = w.wcs_world2pix(rastar, decstar, 0)
		except KeyError: # If there is no WCS solution
			try:
				WCScenter = WCScenter # default to the previous WCS center location
			except NameError: # but, if it's the first image, stop
				print "Initial image has no WCS solution! Stopping."
				sys.exit()

		centerloc = HowellCenter(image, WCScenter, L)
		xcenter[i] = centerloc[0]
		ycenter[i] = centerloc[1]

		photaper = pu.CircularAperture(centerloc, r=rphot)
		bkgaper = pu.CircularAnnulus(centerloc, r_in=rbkgin, r_out=rbkgout)
		raw_table = pu.aperture_photometry(image, photaper)
		bkg_table = pu.aperture_photometry(image, bkgaper)
		bkgmean = bkg_table['aperture_sum'] / bkgaper.area()
		photbkg = bkgmean*photaper.area()
		flux[i] = raw_table['aperture_sum'] - photbkg
		error[i] = numpy.sqrt(raw_table['aperture_sum'] + photbkg)
		bar.update(i+1)
	bar.finish()

	toreturn = numpy.column_stack((jdutc, flux, error, airmass, xcenter, ycenter))
	return toreturn