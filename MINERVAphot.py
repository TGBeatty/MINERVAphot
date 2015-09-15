import numpy, progressbar, astropy, math, sys, os, aplpy, subprocess, json
import photutils as pu
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def MakeBias(firstarg, uselist=False, skip=False, outloc='./', telno='none'):
	if not os.path.exists(outloc):
		os.makedirs(outloc)

	if uselist:
		fitslist = firstarg
		biasloc = outloc+"mbias.fits"
		print "Making bias image..."
		nimages = len(fitslist)
	else:
		night = firstarg
		if telno=='none':
			print "ERROR: You need to specify what telescope in the bias creation. Stopping."
			sys.exit()
		subprocess.call("ls ./"+night+"/*"+telno+".Bias.* > biases.list", shell=True)
		nfiles = len(open("biases.list").readlines())
		if nfiles==0:
			print "ERROR: There are no biases for "+telno+" on the night specified. Stopping."
			sys.exit()
		fitslist = numpy.loadtxt('biases.list', dtype='object')
		biasloc = "./"+night+"/mbias."+telno+".fits"
		print "Making "+telno+" bias image..."
		nimages = len(fitslist)

	if not skip:
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
	else:
		print "Skipping!"

	return biasloc

def MakeDark(firstarg, mbiasloc, uselist=False, skip=False, outloc='./', telno='none'):
	if not os.path.exists(outloc):
		os.makedirs(outloc)

	if uselist:
		fitslist = firstarg
		darkloc = outloc+"mdark.fits"
		print "Making dark image..."
		nimages = len(fitslist)
	else:
		night = firstarg
		if telno=='none':
			print "ERROR: You need to specify what telescope in the dark creation. Stopping."
			sys.exit()
		subprocess.call("ls ./"+night+"/*"+telno+".Dark.* > darks.list", shell=True)
		nfiles = len(open("darks.list").readlines())
		if nfiles==0:
			print "ERROR: There are no darks for "+telno+" on the night specified. Stopping."
			sys.exit()
		fitslist = numpy.loadtxt('darks.list', dtype='object')
		darkloc = "./"+night+"/mdark."+telno+".fits"
		print "Making "+telno+" dark image..."
		nimages = len(fitslist)

	if not skip:
		hdulist = astropy.io.fits.open(mbiasloc)
		bias = hdulist[0].data
		bar = progressbar.ProgressBar(maxval=nimages, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		for i in range(nimages):
			hdulist = astropy.io.fits.open(fitslist[i])
			if i == 0:
				size = (nimages,hdulist[0].data.shape[0],hdulist[0].data.shape[1])
				tomedian = numpy.zeros(size)
			exptime = numpy.float(hdulist[0].header['EXPTIME'])
			tomedian[i,:,:] = (hdulist[0].data - bias) / exptime
			bar.update(i+1)
		bar.finish()

		hdulist[0].data = numpy.median(tomedian, axis=0)
		hdulist.writeto(darkloc, clobber=True)
	else:
		print "Skipping!"

	return darkloc

def MakeFlat(firstarg, mbiasloc, mdarkloc, imgfilter='none', uselist=False, skip=False, outloc='./', telno='none'):
	if not os.path.exists(outloc):
		os.makedirs(outloc)

	if uselist:
		print "Making flat image..."
		if not skip:
			fitslist = firstarg
			flatloc = outloc+"mflat.fits"
			nimages = len(fitslist)
	else:
		night = firstarg
		flatloc = "./"+night+"/mflat."+imgfilter+"."+telno+".fits"
		if telno=='none':
			print "ERROR: You need to specify what telescope to use (T1, T2, etc.) in the flat creation. Stopping."
			sys.exit()
		if imgfilter=='none':
			print "ERROR: You need to specify what filter to use (V, gp, etc.) in the flat creation. Stopping."
			sys.exit()

		if not skip:
			subprocess.call("ls ./"+night+"/*"+telno+".SkyFlat."+imgfilter+".* > flats.list", shell=True)
			nfiles = len(open("flats.list").readlines())
			if nfiles==0:
				AlreadyThere = os.path.isfile(flatloc)
				if AlreadyThere: 
					print "No flats in "+imgfilter+" for "+telno+" on the night of "+night+", but there's already a median flat present. Using that and skipping creation."
					skip = True
				else: 
					print "ERROR: There are no flats in "+imgfilter+" for "+telno+" on the night specified, and no provided median flat. Stopping."
					sys.exit()
			if not skip:
				fitslist = numpy.loadtxt('flats.list', dtype='object')
				nimages = len(fitslist)
		print "Making "+telno+" "+imgfilter+" flat image..."

	if not skip:
		hdulist = astropy.io.fits.open(mbiasloc)
		bias = hdulist[0].data
		hdulist = astropy.io.fits.open(mdarkloc)
		dark = hdulist[0].data
		bar = progressbar.ProgressBar(maxval=nimages, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		for i in range(nimages):
			hdulist = astropy.io.fits.open(fitslist[i])
			if i == 0:
				size = (nimages,hdulist[0].data.shape[0],hdulist[0].data.shape[1])
				tomedian = numpy.zeros(size)
			exptime = numpy.float(hdulist[0].header['EXPTIME'])
			scaleddark = dark * exptime
			tomedian[i,:,:] = (hdulist[0].data - bias - scaleddark) / numpy.median(hdulist[0].data)
			bar.update(i+1)
		bar.finish()

		hdulist[0].data = numpy.median(tomedian, axis=0)
		hdulist.writeto(flatloc, clobber=True)
	else:
		print "Skipping!"

	return flatloc

def BiasDarkFlat(firstarg, mbiasloc, mdarkloc, mflatloc, imgfilter='none', uselist=False, skip=False, outloc='./', telno='none'):
	if not os.path.exists(outloc):
		os.makedirs(outloc)

	if uselist:
		print "Calibrating images..."
		if not skip:
			fitslist = firstarg
			nimages = len(fitslist)
			procimgloc = numpy.empty(nimages, dtype='object') # What will be the list of processed images
	else:
		night = firstarg
		outloc = "./"+night+"/Processed/"
		if not os.path.exists(outloc):
			os.makedirs(outloc)
		if telno=='none':
			print "ERROR: You need to specify what telescope to use (T1, T2, etc.) in the calibration. Stopping."
			sys.exit()
		if imgfilter=='none':
			print "ERROR: You need to specify what filter to use (V, gp, etc.) in the caliberation stage. Stopping."
			sys.exit()

		subprocess.call("find ./"+night+"/ -maxdepth 1 -type f -not -name \"*SkyFlat*\" -name \"*"+telno+".*."+imgfilter+"*\" > images.list", shell=True) # excludes the skyflat images
		nfiles = len(open("images.list").readlines())
		if nfiles==0:
			print "WARNING: There are no science images for a detected filter. Something's wrong somewhere. Continuing."
			skip = True
		fitslist = numpy.loadtxt('images.list', dtype='object')
		nimages = len(fitslist)
		procimgloc = numpy.empty(nimages, dtype='object') # What will be the list of processed images
		print "Calibrating "+telno+" "+imgfilter+" images..."

	if not skip:
		hdulist = astropy.io.fits.open(mbiasloc)
		bias = hdulist[0].data
		hdulist = astropy.io.fits.open(mdarkloc)
		dark = hdulist[0].data
		hdulist = astropy.io.fits.open(mflatloc)
		flat = hdulist[0].data 

		bar = progressbar.ProgressBar(maxval=nimages, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		for i in range(nimages):
			hdulist = astropy.io.fits.open(fitslist[i])
			image = hdulist[0].data
			exptime = numpy.float(hdulist[0].header['EXPTIME'])
			scaleddark = dark * exptime
			hdulist[0].data = (image - bias - scaleddark) / flat
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

def CalibrateNight(night, telno, skipbias=False, skipdark=False, skipflat=False, skipproc=False):
	print "-----------------------------------------------------------------------"
	print "Calibrating an entire night at once, based on the target file."
	print "Night: "+night
	print "Telescope: "+telno
	print "-----------------------------------------------------------------------"

	targetfilename = night+'.txt'

	readin = []
	with open(targetfilename, 'r') as targetfile:
		for line in targetfile:
			readin.append(json.loads(line))
	startcals = readin[0]
	endcals = readin[1]
	ntargetsOG = len(readin)-2
	nameOG = numpy.zeros(ntargetsOG, dtype='object')
	raOG = numpy.zeros(ntargetsOG, dtype='object')
	decOG = numpy.zeros(ntargetsOG, dtype='object')
	exptimeOG = numpy.zeros(ntargetsOG, dtype='object')
	filtersOG = numpy.zeros(ntargetsOG, dtype='object')
	for i in range(ntargetsOG):
		targline = readin[i+2]
		nameOG[i] = targline['name']
		raOG[i] = targline['ra']
		decOG[i] = targline['dec']
		exptimeOG[i] = targline['exptime']
		filtersOG[i] = targline['filter']

	# find the unique filters
	filtersIndiv = numpy.zeros(1, dtype='object')
	for i in range(ntargetsOG):
		nfilter = len(filtersOG[i])
		filtertemp = filtersOG[i]
		for j in range(nfilter): filtersIndiv = numpy.append(filtersIndiv, filtertemp[j])
	filtersIndiv = numpy.delete(filtersIndiv, 0)
	filters = numpy.unique(filtersIndiv)

	biasloc = MakeBias(night, telno=telno, skip=skipbias)
	darkloc = MakeDark(night, biasloc, telno=telno, skip=skipdark)
	procimgs = {}
	for imgfilter in filters:
		flatloc = MakeFlat(night, biasloc, darkloc, imgfilter=imgfilter, telno=telno, skip=skipflat)
		procimgs[imgfilter] = BiasDarkFlat(night, biasloc, darkloc, flatloc, imgfilter=imgfilter, telno=telno, skip=skipproc)

	return procimgs, filters

def CalibrateList(unproclist, biaslist, darklist, flatlist, outloc, skipbias=False, skipdark=False, skipflat=False, skipproc=False):
	print "-----------------------------------------------------------------------"
	print "Calibrating a list of images, using lists of calibration images."
	print "-----------------------------------------------------------------------"

	biasloc = MakeBias(biaslist, outloc=outloc, uselist=True, skip=skipbias)
	darkloc = MakeDark(darklist, biasloc, outloc=outloc, uselist=True, skip=skipdark)
	flatloc = MakeFlat(flatlist, biasloc, darkloc, outloc=outloc, uselist=True, skip=skipflat)
	procimgs = {}
	procimgs['nofilter'] = BiasDarkFlat(unproclist, biasloc, darkloc, flatloc, outloc=outloc, uselist=True, skip=skipproc)

	return procimgs


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

	toreturn = numpy.column_stack((jdutc, airmass, flux, error, xcenter, ycenter))
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

	toreturn = numpy.column_stack((jdutc, airmass, flux, error, xcenter, ycenter))
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
		print "The first aperture will be the target."
		print "Keystrokes need to happen with the plot window active."
		print "Left-click: place aperture center"
		print "Backspace: delete previous aperture"
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

def SimpleDetrend(photresults):
	ncols = photresults.shape[1]
	nstars = (ncols-2)/4

	data = {}
	data['time'] = photresults[:,0]
	data['airmass'] = photresults[:,1]
	errstart = 2+nstars
	xcstart = errstart + nstars
	ycstart = xcstart + nstars
	for i in range(nstars):
		data["flux"+str(i)] = photresults[:,(2+i)]
		data["error"+str(i)] = photresults[:,(errstart+i)]
		data["xcent"+str(i)] = photresults[:,(xcstart+i)]
		data["ycent"+str(i)] = photresults[:,(ycstart+i)]

	target = [0]
	detrend = numpy.arange(1,nstars)

	#print "Detrend stars are star numbers:"
	#print detrend

	targetflux = photresults[:,(2+target[0])]
	detrendflux = photresults[:,(2+detrend)]
	targeterror = photresults[:,(errstart+target[0])]
	detrenderror = photresults[:,(errstart+detrend)]

	normtarget = targetflux / numpy.median(targetflux)
	normtargeterror = targeterror / numpy.median(targetflux)
	detrendsum = numpy.sum(detrendflux, axis=1)
	detrendsumerror = numpy.sqrt(numpy.sum(detrenderror**2., axis=1))
	normdetrend = detrendsum / numpy.median(detrendsum)
	normdetrenderror = detrendsumerror / numpy.median(detrendsum)

	targetdetrend = normtarget / normdetrend
	targetdetrenderror = targetdetrend * numpy.sqrt((normtargeterror/normtarget)**2. + (normdetrenderror/normdetrend)**2.)

	return targetdetrend, targetdetrenderror
