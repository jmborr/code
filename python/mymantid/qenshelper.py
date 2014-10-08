import numpy

def normalizeHistogram(inputworkspace, outputworkspace=None, RangeLower=None, RangeUpper=None):
	'''Set to unity the integral of the histogram. NOTE: the bin width is used as dx in the integral'''
	kwargs = {}
	if RangeLower: kwargs['RangeLower']=RangeLower
	if RangeUpper: kwargs['RangeUpper']=RangeUpper
	if not outputworkspace: outputworkspace = inputworkspace
	ConvertToHistogram(InputWorkspace=inputworkspace,OutputWorkspace='junk_normalizeHistogram')
	Integration(InputWorkspace='junk_normalizeHistogram',OutputWorkspace='junk2_normalizeHistogram', **kwargs)
	# Recall that Integration algorithm "fails" to multiply by the bin width
	ws=mtd['junk_normalizeHistogram']
	dX = ws.dataX(0)[1] - ws.dataX(0)[0]
	Scale(InputWorkspace='junk2_normalizeHistogram', Factor=dX,Operation='Multiply',OutputWorkspace='junk2_normalizeHistogram')
	Divide(LHSWorkspace=inputworkspace,RHSWorkspace='junk2_normalizeHistogram',OutputWorkspace=outputworkspace)
	DeleteWorkspace(Workspace='junk_normalizeHistogram')
	DeleteWorkspace(Workspace='junk2_normalizeHistogram')

def setToUnity(iwsg):
	'''Set I(Q,t=0)=1'''
	ExtractSingleSpectrum(InputWorkspace='{0}_fq0'.format(iwsg),WorkspaceIndex=0,OutputWorkspace='junk')
	Transpose(InputWorkspace='junk', OutputWorkspace='junk')
	Divide(LHSWorkspace='{0}_fqt.Re'.format(iwsg), RHSWorkspace='junk', OutputWorkspace='{0}_fqt.Re'.format(iwsg))
	DeleteWorkspace(Workspace='junk')

def removeFlatBackground(sqw, outputworkspace=None):
	'''S(Q,E) goes to a constant as E goes to infinity, instead of zero'''
	import numpy
	if not outputworkspace:
		outputworkspace = sqw
	else:
		CloneWorkspace(InputWorkspace=sqw, OutputWorkspace=outputworkspace)
	percent=0.05
	ws=mtd[outputworkspace]
	for ix in range( ws.getNumberHistograms() ):
		print "len=",len(ws.dataY(ix))
		N=len(ws.dataY(ix))*percent
		print 'N=',N
		average = numpy.sum(ws.dataY(ix)[-N:])/N
		print 'average=',average
		ws.dataY(ix)[:] -= average

def rebinInQ(wsname, params):
	'''Rebin along the vertical axes'''
	Transpose(InputWorkspace=wsname, OutputWorkspace=wsname)
	Rebin(InputWorkspace=wsname, Params=params, OutputWorkspace=wsname)
	Transpose(InputWorkspace=wsname, OutputWorkspace=wsname)

def clipIQ0(wsname, outname=None, mode='quadratic', **kwargs):
	'''Remove the peak of IQt. We don't know the shape of I(Q,t) in the range [0,dt], 
	where dt is the time between simulation frames. This produces a high frequency background in S(Q,E).
	The background can be removed if the extrapolate I(Q,t=0) from a fit to the nearby points I(Q,t>0).
	Since I(Q,t=0) must be one, one can always renormalize the intensities so that I(Q,0)=1, although
	this is not done in this routine
	'''
	if outname:
		CloneWorkspace(InputWorkspace=wsname, OutputWorkspace=outname)
	else:
		outname=wsname
	ws = mtd[outname] # print 'wsname=', wsname
	nk = ws.getNumberHistograms() #print 'nk=',nk
	if mode=='quadratic':
		'''fit I(Q,0) with quadratic function'''
		nt=4
		if 'nt' in kwargs.keys(): nt=kwargs['nt']
		for k in range(nk):
			Y = ws.dataY(k)
			hY = numpy.argmax(Y) # index with maximum of I(Q,t)
			X = ws.dataX(k)
			x = X[hY+1:hY+1+nt]
			if len(X) != len(Y):
				#Histogram data, the maximum may be duplicated if I(Q,t) is even. Find the highest index then
				hY = numpy.where( Y==Y[hY] )[-1]
				x = ( X[hY+1:hY+1+nt] + X[hY+2:hY+2+nt] ) / 2.0
			y = Y[hY+1:hY+1+nt]
			p = numpy.polyfit( x, y, 2)
			Y[ numpy.where( Y==Y[hY] ) ] = p[-1] # substitute the maximum intensity with the intercept of the quadratic fit


def clipSQ0(wsname, outname=None, mode='truncate', **kwargs):
	'''Remove the peak of SQE
	mode
		truncate: clipp the peak, leaving a flat central area
		lorentzian: fit a lorentzian and substite the peak for the lorentzian values
	'''
	import numpy
	if outname:
		CloneWorkspace(InputWorkspace=wsname, OutputWorkspace=outname)
	else:
		outname=wsname
	ws = mtd[outname] # print 'wsname=', wsname
	nk = ws.getNumberHistograms() #print 'nk=',nk
	if mode=='truncate':
		ratio = 0.1
		if 'ratio' in kwargs.keys(): ratio = kwargs['ratio']
		for k in range(nk):
			Y = ws.dataY(k)
			nx = len(Y) # print 'nx=',nx
			x0 = numpy.argmax(Y)  # index of histogram maximum print 'x0=',x0
			x0l = [x0,] # list of index where to clip SQE
			y0 = Y[x0] # value of histogram maximum print 'y0=', y0
			# Search adjacent points x < x0
			search = 1
			y0left = y0
			x = x0 - 1
			while search and x >= 0:
				if Y[x] > y0*ratio:
					x0l.append(x)
					y0left = Y[x]
					x = x - 1
				else:
					yleft = Y[x]
					search = 0
			search = 1
			y0right = y0
			x = x0 + 1
			while search and x <= nx:
				if Y[x] > y0*ratio:
					x0l.append(x)
					y0right = Y[x]
					x = x + 1
				else:
					yright = Y[x]
					search = 0
			# clip values
			yfinal = (yleft+yright) / 2 # print 'yfinal=', yfinal  print 'x0l=',str(x0l)
			for x in x0l:
				Y[x] = yfinal
	elif mode == 'lorentzian':
		def lorentz(x, intensity, FWHM):
			''' Lorentzian function '''
			return intensity / (FWHM**2 + x*x)
		from scipy.optimize import curve_fit
		emin = 0.0034 # assumed Basis resolution
		emax = 5*emin
		if 'emin' in kwargs.keys(): emin = kwargs['emin']
		if 'emax' in kwargs.keys(): emax = kwargs['emax']
		for k in range(nk):
			Y = ws.dataY(k)
			nY = len(Y)
			X = ws.dataX(k)
			nX = len(X)
			indexes = numpy.intersect1d( numpy.where(X>emin)[0], numpy.where(X<emax)[0] )
			y = Y[indexes]
			x = X[indexes]
			if nX != nY:
				x = (X[indexes] + X[indexes+1])/2  # deal with histogram data
			popt, pcov = curve_fit(lorentz, x, y, p0=[y[0],emin]) # Fit (x,y) with a Lorentzian1D
			# Evaluate the lorentzian in [-emin, emin] and substitute the data with it
			x0 = X[indexes[0]+1]  #extend the region one bin beyond the emin threshold
			indexes = numpy.intersect1d( numpy.where(X>-x0)[0], numpy.where(X<x0)[0] )
			x = X[indexes]
			if nX != nY:
				x = (X[indexes] + X[indexes+1])/2  # deal with histogram data
			y = lorentz(x, *popt)
			Y[indexes] = y	

def cancel_negative_intensities(inputworkspace, outputworkspace=None):
	'''Navigate the structure factor and cancel negative intensities by taking positive intensities from other energy values.
	This transformation preserves the integrated intensity'''
	if outputworkspace: 
		CloneWorkspace(InputWorkspace=inputworkspace, OutputWorkspace=outputworkspace)
	else:
		outputworkspace=inputworkspace
	ws = mtd[outputworkspace] # print 'wsname=', wsname
	nk = ws.getNumberHistograms() #print 'nk=',nk
	for k in range(nk):
		print 'k=',k
		Y = ws.dataY(k)
		en = len(Y) # number of energy points
		eh = numpy.argmax(Y) # index where maximum of SQE is located
		#Smooth the negative energies
		area=0.0 #intensity "to the left" of ie
		for ie in range(0,eh):
                        #print 'k=',k,'area=',area,'ie=',ie
			# Negative intensity !
			if Y[ie] < 0.0:
                                #print 'ie',ie,'Y[ie]=',Y[ie]
				# smooth Y[ie] grabbing intensity from Y[je<ie]	
				if area > 0.0:
					for je in range(0,ie):
						if Y[je] > 0.0:
							# intensity at Y[je] not enough to set Y[ie] to zero
							if Y[je] < abs(Y[ie]):
								Y[ie] += Y[je]
								area -= Y[je]
								Y[je] = 0.0
							else:
								Y[je] += Y[ie]
								area  += Y[ie]
								Y[ie] = 0.0
						if Y[ie] == 0.0:
							break
			# If there was not enough area in Y[je<ie], grab intensity from Y[ie<je]
			if Y[ie] < 0.0:
				for je in range(ie+1,eh):
					if Y[je] > 0.0:
						# intensity at Y[je] not enough to set Y[ie] to zero
						if Y[je] < abs(Y[ie]):
							Y[ie] += Y[je]
							Y[je] = 0.0
						else:
							Y[je] += Y[ie]
							Y[ie] = 0.0
					if Y[ie] == 0.0:
						break
			else:
				area += Y[ie] # Y[ie] could have been positive from the start
		#Smooth the positive energies
		area=0.0 #intensity "to the right" of ie
		for ie in range(en-1, eh, -1):
                        #print 'k=',k,'area=',area,'ie=',ie
			# Negative intensity !
			if Y[ie] < 0.0:
				# smooth Y[ie] grabbing intensity from Y[ie<je]	
				if area > 0.0:
					for je in range(en-1,ie,-1):
						if Y[je] > 0.0:
							# intensity at Y[je] not enough to set Y[ie] to zero
							if Y[je] < abs(Y[ie]):
								Y[ie] += Y[je]
								area -= Y[je]
								Y[je] = 0.0
							else:
								Y[je] += Y[ie]
								area  += Y[ie]
								Y[ie] = 0.0
						if Y[ie] == 0.0:
							break
			# If there was not enough area in Y[je<ie], grab intensity from Y[je<ie]
			if Y[ie] < 0.0:
				for je in range(ie-1,eh,-1):
					if Y[je] > 0.0:
						# intensity at Y[je] not enough to set Y[ie] to zero
						if Y[je] < abs(Y[ie]):
							Y[ie] += Y[je]
							Y[je] = 0.0
						else:
							Y[je] += Y[ie]
							Y[ie] = 0.0
					if Y[ie] == 0.0:
						break
			else:
				area += Y[ie] # Y[ie] could have been positive from the start

