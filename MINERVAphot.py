import numpy, progressbar, astropy, math, sys, os, aplpy
import photutils as pu
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def MakeBias(fitslist, outloc, skip=False):
	if not os.path.exists(outloc):
		os.makedirs(outloc)

	nimages = len(fitslist)
	biasloc = outloc+"mbias.fits"

	if not skip:
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
		hdulist.writeto(biasloc, clobber=True)

	return biasloc

def MakeDark(fitslist, outloc, mbiasloc, skip=False):
	if not os.path.exists(outloc):
		os.makedirs(outloc)

	nimages = len(fitslist)
	darkloc = outloc+"mdark.fits"

	if not skip:
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
		hdulist.writeto(darkloc, clobber=True)

	return darkloc

def MakeFlat(fitslist, outloc, mbiasloc, mdarkloc, skip=False):
	if not os.path.exists(outloc):
		os.makedirs(outloc)

	nimages = len(fitslist)
	flatloc = outloc+"mflat.fits"

	if not skip:
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
		hdulist.writeto(flatloc, clobber=True)

	return flatloc

def BiasDarkFlat(fitslist, outloc, mbiasloc, mdarkloc, mflatloc, skip=False):
	if not os.path.exists(outloc):
		os.makedirs(outloc)

	nimages = len(fitslist)
	procimgloc = numpy.empty(nimages, dtype='object') # What will be the list of processed images

	if not skip:
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
	else:
		for i in range(nimages):
			splitname = fitslist[i].split('/')
			filenamesplit = splitname[-1].partition('.fit')
			procimgloc[i] = outloc+filenamesplit[0]+".proc"+filenamesplit[1]+filenamesplit[2]

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

def SingleAper(fitslist, AperProps):

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
			WCScenter = w.wcs_world2pix(AperProps['rastar']*15., AperProps['decstar'], 0)
		except KeyError: # If there is no WCS solution
			try:
				WCScenter = WCScenter # default to the previous WCS center location
			except NameError: # but, if it's the first image, stop
				print "Initial image has no WCS solution! Stopping."
				sys.exit()
		if (WCScenter[0]<0.) or (WCScenter[0]>2048.): 
			print "The WCS solution has given an invalid pixel location! Stopping."
			print WCScenter
			print "Image name is:", fitslist[i]
			sys.exit()
		if (WCScenter[1]<0.) or (WCScenter[1]>2048.): 
			print "The WCS solution has given an invalid pixel location! Stopping."
			print WCScenter
			print "Image name is:", fitslist[i]
			sys.exit()

		centerloc = HowellCenter(image, WCScenter, AperProps['L'])
		xcenter[i] = centerloc[0]
		ycenter[i] = centerloc[1]

		photaper = pu.CircularAperture(centerloc, r=AperProps['rAp'])
		bkgaper = pu.CircularAnnulus(centerloc, r_in=AperProps['rBkgIn'], r_out=AperProps['rBkgOut'])
		raw_table = pu.aperture_photometry(image, photaper)
		bkg_table = pu.aperture_photometry(image, bkgaper)
		bkgmean = bkg_table['aperture_sum'] / bkgaper.area()
		photbkg = bkgmean*photaper.area()
		flux[i] = raw_table['aperture_sum'] - photbkg
		error[i] = numpy.sqrt(raw_table['aperture_sum'] + photbkg)
		bar.update(i+1)
		hdulist.close()
	bar.finish()

	toreturn = numpy.column_stack((jdutc, flux, error, airmass, xcenter, ycenter))
	return toreturn

def MultiAper(fitslist, AperProps):

	nstars = AperProps['rastar'].size
	RAlist = AperProps['rastar']
	Declist = AperProps['decstar']

	nimages = len(fitslist)
	jdutc = numpy.zeros(nimages)
	flux = numpy.zeros((nimages,nstars,))
	error = numpy.zeros((nimages,nstars))
	xcenter = numpy.zeros((nimages,nstars))
	ycenter = numpy.zeros((nimages,nstars))
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

		w = astropy.wcs.WCS(hdulist[0].header)
		for j in range(nstars):
			try:
				wcscheck = hdulist[0].header['WCSAXES']
				WCScenter = w.wcs_world2pix(RAlist[j]*15., Declist[j], 0)
			except KeyError: # If there is no WCS solution
				try:
					WCScenter = WCScenter # default to the previous WCS center location
				except NameError: # but, if it's the first image, stop
					print "Initial image has no WCS solution! Stopping."
					sys.exit()
			if (WCScenter[0]<0.) or (WCScenter[0]>2048.): 
				print "The WCS solution has given an invalid pixel location! Stopping."
				print WCScenter
				print "Image name is:", fitslist[i]
				sys.exit()
			if (WCScenter[1]<0.) or (WCScenter[1]>2048.): 
				print "The WCS solution has given an invalid pixel location! Stopping."
				print WCScenter
				print "Image name is:", fitslist[i]
				sys.exit()

			centerloc = HowellCenter(image, WCScenter, AperProps['L'])
			xcenter[i,j] = centerloc[0]
			ycenter[i,j] = centerloc[1]

			photaper = pu.CircularAperture(centerloc, r=AperProps['rAp'])
			bkgaper = pu.CircularAnnulus(centerloc, r_in=AperProps['rBkgIn'], r_out=AperProps['rBkgOut'])
			raw_table = pu.aperture_photometry(image, photaper)
			bkg_table = pu.aperture_photometry(image, bkgaper)
			bkgmean = bkg_table['aperture_sum'] / bkgaper.area()
			photbkg = bkgmean*photaper.area()
			flux[i,j] = raw_table['aperture_sum'] - photbkg
			error[i,j] = numpy.sqrt(raw_table['aperture_sum'] + photbkg)
		hdulist.close()
		bar.update(i+1)
	bar.finish()

	toreturn = numpy.column_stack((jdutc, flux, error, airmass, xcenter, ycenter))
	return toreturn

def on_click(event):
	global PickedLocsPix, w, fig, WCSLocsStore, FirstFlag, image

	if event.inaxes is not None:
		recentered = HowellCenter(image, [event.xdata,event.ydata], 20)
		if FirstFlag: 
			WCSLocsStore = w.wcs_pix2world(recentered[0], recentered[1], 0)
			PickedLocsPix = recentered
			fig.show_markers(WCSLocsStore[0], WCSLocsStore[1], layer='Markers')
			FirstFlag = False
		else:
			WCSLocsStore = numpy.vstack((WCSLocsStore,w.wcs_pix2world(recentered[0], recentered[1], 0)))
			PickedLocsPix = numpy.vstack((PickedLocsPix,recentered))
			fig.show_markers(WCSLocsStore[:,0], WCSLocsStore[:,1], layer='Markers')
		fig.refresh()
	else:
		print 'Error: clicked outside axes bounds but inside plot window'

def on_key(event):
	global PickedLocsPix, fig, WCSLocsStore, FirstFlag

	if event.key=='backspace':
		if not FirstFlag:
			if WCSLocsStore.size > 0:
				WCSLocsStore = numpy.delete(WCSLocsStore, -1, 0)
				PickedLocsPix = numpy.delete(PickedLocsPix, -1, 0)
				if WCSLocsStore.size > 0:
					fig.show_markers(WCSLocsStore[:,0], WCSLocsStore[:,1], layer='Markers')
				else: fig.remove_layer('Markers') # show nothing if no entries
			fig.refresh()
	if event.key=='enter':
		plt.close()
		plt.clf()
	if event.key=='escape':
		plt.close()
		plt.clf()

def GraphicalPickLocs(imagename):
	global PickedLocsPix, w, fig, WCSLocsStore, FirstFlag, image

	PickedLocsPix = numpy.zeros((2,1))
	WCSLocsStore = numpy.zeros((2,1))
	FirstFlag = True
	StillNeedToPick = True

	hdulist = astropy.io.fits.open(imagename)
	image = hdulist[0].data 
	w = astropy.wcs.WCS(hdulist[0].header)
	while StillNeedToPick:
		mplfig = plt.figure()
		fig = aplpy.FITSFigure(imagename, figure=mplfig)
		fig.show_grayscale()
		fig.set_xaxis_coord_type('longitude')
		fig.set_yaxis_coord_type('latitude')
		fig.tick_labels.set_xformat('dd.dddddd')
		fig.tick_labels.set_yformat('dd.dddddd')
		cid1 = mplfig.canvas.callbacks.connect('button_press_event', on_click)
		cid2 = mplfig.canvas.callbacks.connect('key_press_event', on_key)
		print "-----------------------------------------------------------------------"
		print "Pick locations for the apertures. Don't worry, they will be recentered!"
		print "Keystrokes need to happen with the plot window active."
		print "Left-click: place aperture center"
		print "Enter / escape: finish and close window"
		print "-----------------------------------------------------------------------"
		plt.show()
		plt.clf()
		plt.close()
		mplfig.canvas.mpl_disconnect(cid1)
		mplfig.canvas.mpl_disconnect(cid2)
		if (FirstFlag or len(PickedLocsPix)==0): # didn't pick at all
			print "ERROR: you didn't pick any locations. Trying again..."
		else:
			try: 
				test = PickedLocsPix.shape
				PickedLocsWCS = w.wcs_pix2world(PickedLocsPix[:,0], PickedLocsPix[:,1]	, 0)
				StillNeedToPick = False
			except AttributeError: # picked just one location
				PickedLocsWCS = w.wcs_pix2world(PickedLocsPix[0], PickedLocsPix[1], 0)
				StillNeedToPick = False

	hdulist.close()
	return PickedLocsWCS