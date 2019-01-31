import os.path
import numpy as np
import utils
import processor 
import scipy.integrate as intg

GRIBflag = True
MERRAflag = True

try:
    import era
except:
    GRIBflag = False

try:
    import narr
except:
    GRIBflag = False

try:
    import merra
except:
    MERRAflag = False

import scipy.interpolate as si
import matplotlib.pyplot as plt
import sys

def sanity_check(model):
    if model in ('ECMWF','ERA','NARR') and GRIBflag==False:
        print 'No module name pygrib found.'
        print 'To use ECMWF and NARR, please install pygrib.'

    if model=='MERRA' and MERRAflag==False:
        print 'No module name pyhdf found'
        print 'To use MERRA, please install pyhdf'
        print 'Visit: http://pysclint.sourceforge.net/pyhdf'


##############Creating a class object for PyAPS use.
class PyAPS_rdr:
	'''Class for dealing with Atmospheric phase corrections in radar coordinates.
	Operates on one weather model file and one Geo.rsc file at a time.'''
	
	def __init__(self,gribfile,demfile,grib='ECMWF',humidity='Q',demtype=np.float32,demfmt='RMG',verb=False, Del='comb'):
		'''Initiates the data structure for atmos corrections in geocoded domain.'''
                sanity_check(grib)
                grib = grib.upper()
		humidity = humidity.upper()
		demfmt = demfmt.upper()
		
		assert grib in ('ERA','NARR','ECMWF','MERRA'), 'PyAPS: Undefined grib file source.'
		self.grib = grib
		
		assert humidity in ('Q','R'), 'PyAPS: Undefined humidity.'
		self.humidity = humidity
		
		if self.grib in ('NARR','MERRA'):
				assert self.humidity in ('Q'), 'PyAPS: Relative humidity not provided by NARR/MERRA.'

		assert os.path.isfile(gribfile), 'PyAPS: GRIB File does not exist.'
		self.gfile = gribfile
		
		assert os.path.isfile(demfile), 'PyAPS: DEM file does not exist.'
		self.hfile = demfile

		assert demfmt in ('RMG','HGT'), 'PyAPS: DEM Format can be RMG or HGT'
		self.fmt = demfmt
	
		self.demtype = demtype
			
		self.dict = processor.initconst()

		[lon,lat,nx,ny,dpix] = utils.rd_rsc(self.hfile,verbose=verb)

		if grib in ('ERA','ECMWF'):
			if GRIBflag:
				self.hgtscale=((self.dict['maxAlt']-self.dict['minAlt'])/self.dict['nhgt'])/0.703
				self.rdrscale = 70000./dpix
				self.bufspc = 1.2
			else:
				print '================================'
                                print '********************************'
                                print '  pyGrib needs to be installed  '
                                print '    No ECMWF or NARR possible   '
                                print '********************************'
				print '================================'
				sys.exit(1)
		elif grib in ('NARR'):
			if GRIBflag:
                                self.hgtscale=((self.dict['maxAlt']-self.dict['minAlt'])/self.dict['nhgt'])/0.703
                                self.rdrscale = 70000./dpix
                                self.bufspc = 1.2
                        else:
                                print '================================'
                                print '********************************'
                                print '  pyGrib needs to be installed  '
                                print '    No ECMWF or NARR possible   '
                                print '********************************'
                                print '================================'
                                sys.exit(1)
			self.hgtscale = ((self.dict['maxAlt']-self.dict['minAlt'])/self.dict['nhgt'])/0.3
			self.rdrscale = 32000./dpix
			self.bufspc = 1.2
		elif grib in ('MERRA'):
			if MERRAflag:
				self.hgtscale = ((self.dict['maxAlt']-self.dict['minAlt'])/self.dict['nhgt'])/0.5
				self.rdrscale = 32000./dpix
				self.bufspc = 1.0 
			else:
				print '================================'
				print '********************************'
				print '  pyHdf needs to be installed   '
				print '        No MERRA possible       '
				print '********************************'
				print '================================'
				sys.exit(1)
			
		######Problems near the international dateline
		if grib in ('MERRA'):
			self.minlon = lon.min()-self.bufspc
			self.maxlon = lon.max()+self.bufspc
			self.minlat = lat.min()-self.bufspc
			self.maxlat = lat.max()+self.bufspc 
			self.nx = nx
			self.ny = ny
		else:	
			lon[lon < 0.] += 360.0
			self.minlon = lon.min()-self.bufspc
			self.maxlon = lon.max()+self.bufspc
			self.minlat = lat.min()-self.bufspc
			self.maxlat = lat.max()+self.bufspc
			self.nx = nx
			self.ny = ny
		
		if self.grib in ('ERA'):
			[lvls,latlist,lonlist,gph,tmp,vpr] = era.get_era(self.gfile,self.minlat,self.maxlat,self.minlon,self.maxlon,self.dict,humidity=self.humidity,verbose=verb)

		elif self.grib in ('ECMWF'):
			[lvls,latlist,lonlist,gph,tmp,vpr] = era.get_ecmwf(self.gfile,self.minlat,self.maxlat,self.minlon,self.maxlon,self.dict,humidity=self.humidity,verbose=verb)
			
		elif self.grib in ('NARR'):
			[lvls,latlist,lonlist,gph,tmp,vpr] = narr.get_narr(self.gfile,self.minlat,self.maxlat,self.minlon,self.maxlon,self.dict,verbose=verb)
		elif self.grib in ('MERRA'):
			[lvls,latlist,lonlist,gph,tmp,vpr] = merra.get_merra(self.gfile,self.minlat,self.maxlat,self.minlon,self.maxlon,self.dict,verbose=verb)
			lonlist[lonlist < 0.] += 360.0 

		[xi,yi] = utils.glob2rdr(nx,ny,lat,lon,latlist,lonlist)

		hgt = np.linspace(self.dict['minAlt'],self.dict['maxAlt'],self.dict['nhgt'])
		
		[Pi,Ti,Vi] = processor.intP2H(lvls,hgt,gph,tmp,vpr,self.dict,verbose=verb)
		
		[DDry,DWet] = processor.PTV2del(Pi,Ti,Vi,hgt,self.dict,verbose=verb)
                if Del in ('comb','Comb'):
                    Delfn = DDry+DWet
                elif Del in ('dry','Dry'):
                    Delfn = DDry
                elif Del in ('wet','Wet'):
                    Delfn = DWet
                else:
                    print 'Unrecognized delay type'
                    sys.exit(1)

		if self.grib in ('NARR'):
			[Delfn,latlist,lonlist] = narr.intdel(hgt,latlist,lonlist,Delfn)

		self.Delfn = Delfn
		self.latlist = latlist
		self.lonlist = lonlist
		self.xi = xi
		self.yi = yi
		self.lat = lat
		self.lon = lon
		self.hgt = hgt
		self.Pi = Pi
		self.Ti = Ti
		self.Vi = Vi
		self.verb = verb

	def getlindelay(self,dataobj,inc=0.0,wvl=4*np.pi):
		'''Write delay to a matrix / HDF5 object or a file directly. LinearNDInterpolator is used, and is not really stable.

		Args:
			* dataobj  (str or HDF5 or np.array): Final output. If str, output is written to file.

		Kwargs:
			* inc  (np.float): Incidence angle in degrees. Default is vertical.
			* wvl  (np.float): Wavelength in meters. Default output results in delay in meters.

		.. note::
			If dataobj is string, output is written to the file.
			If np.array or HDF5 object, it should be of size (ny,nx).'''


		self.rdrfnc = processor.make3dintp(self.Delfn,self.lonlist,self.latlist,self.hgt,self.hgtscale)

		minAltp = self.dict['minAltP']
		xarr = np.arange(1.,self.nx+1.)
		fin = open(self.hfile,'rb');

		outFile = isinstance(dataobj,str)
		
		if outFile:
			fout = open(dataobj,'wb')
		else:
			assert ((dataobj.shape[0]==self.ny) & (dataobj.shape[1]==self.nx)), 'PyAPS: Not a valid data object.' 
			
		cinc = np.cos(inc*np.pi/180.0)
		toto = utils.ProgressBar(maxValue=self.ny)
		for m in range(self.ny):
			if self.fmt in ('HGT'):
				dem = np.fromfile(file=fin,dtype=self.demtype,count=self.nx)
			elif self.fmt in ('RMG'):
				dem = np.fromfile(file=fin,dtype=self.demtype,count=2*self.nx)
				dem = dem[self.nx:]

			dem[dem<minAltp] = minAltp
			demy = dem.astype(np.float64)
			llh = np.zeros((self.nx,3))
			yarr = (m+1)*np.ones((xarr.shape))
			[xin,yin] = utils.rdr2glob(self.nx,self.ny,self.lat,self.lon,xarr,yarr)
			llh[:,0] = xin
			llh[:,1] = yin
			llh[:,2] = demy/self.hgtscale
			res = self.rdrfnc(llh)
			res = res*np.pi*4.0/(cinc*wvl)
			res = res.flatten()
			if outFile:
				resy = res.astype(np.float32)
				resy.tofile(fout)
			else:
				dataobj[m,:] = res
			toto.update(m,every=5)		

		toto.close()

		if outFile:
			fout.close()

		fin.close()

	def getlingeodelay(self,dataobj,lat=None,lon=None,inc=None,wvl=np.pi*4.):
		'''Write delay to a matrix / HDF5 object or a file directly. This is used when the latitude and longitude values are available for each radar pixel. Incidence angle can be a constant or a file name.
		LinearNDinterpolator is used and is not that stable...

		Args:
			* dataobj  (str or HDF5 or np.array): Final output. If str, output is written to file.
		Kwargs:
			* lat  (str)            : Path to the latitude file (np.float32).
			* lon  (str)            : Path to the longitude file (np.float32).
			* inc  (str or np.float): Path to incidence angle file in degrees (str) or a constant float.
			* wvl  (np.float)       : Wavelength in meters. Default output results in delay in meters.

		.. note::
			If dataobj is string, output is written to the file.
			If np.array or HDF5 object, it should be of size (ny,nx).'''

		self.fnc = processor.make3dintp(self.Delfn,self.lonlist,self.latlist,self.hgt,self.hgtscale)

		minAltp = self.dict['minAltP']

		outFile = isinstance(dataobj,str)
		
		if outFile:
			fout = open(dataobj,'wb')
		else:
			assert ((dataobj.shape[0]==self.ny) & (dataobj.shape[1]==self.nx)), 'PyAPS: not a valid data object.'

		assert lat is not None, 'PyAPS: Need a valid latitude file.'
		assert lon is not None, 'PyAPS: Need a valid longitude file.'
	
		if isinstance(inc,float) or isinstance(inc,np.float64) or isinstance(inc,np.float32):
			cinc = np.cos(inc*np.pi/180.)
			incFileFlag = False
		else:
			assert inc is not None, 'PyAPS: Need a valid incidence angle file or constant.'
			incFileFlag = True
			incin = open(inc,'rb')

		latin = open(lat,'rb')
		lonin = open(lon,'rb')
		demin = open(self.hfile,'rb')

		for m in range(self.ny):
			if self.fmt in ('HGT'):
				dem = np.fromfile(file=demin,dtype=self.demtype,count=self.nx)
			elif self.fmt in ('RMG'):
				dem = np.fromfile(file=demin,dtype=self.demtype,count=2*self.nx)
				dem = dem[self.nx:]

			laty = np.fromfile(file=latin,dtype=np.float32,count=self.nx)
			lonx = np.fromfile(file=lonin,dtype=np.float32,count=self.nx)

			lonx[lonx<0] += 360.
			
			if incFileFlag:
				incz = np.fromfile(file=incin,dtype=np.float32,count=self.nx)
				cinc = np.cos(incz*np.pi/180.)

			dem[dem<minAltp] = minAltp
			demy = dem.astype(np.float64)
			llh = np.zeros((self.nx,3))
			llh[:,0] = lonx
			llh[:,1] = laty
			llh[:,2] = demy/self.hgtscale
			res = self.fnc(llh)
			res = res*np.pi*4.0/(cinc*wvl)
			res = res.flatten()
			if outFile:
				resy = res.astype(np.float32)
				resy.tofile(fout)
			else:
				dataobj[m,:] = res

		latin.close()
		lonin.close()
		demin.close()
		if incFileFlag:
			incin.close() 


	def getdelay(self,dataobj,inc=0.0,wvl=4*np.pi):
		'''Write delay to a matrix / HDF5 object or a file directly. Bilinear interpolation at each elevation level is used.

                Args:
                        * dataobj  (str or HDF5 or np.array): Final output. If str, output is written to file.

                Kwargs:
                        * inc  (np.float): Incidence angle in degrees. Default is vertical.
                        * wvl  (np.float): Wavelength in meters. Default output results in delay in meters.

                .. note::
                        If dataobj is string, output is written to the file.
                        If np.array or HDF5 object, it should be of size (ny,nx).'''

		minAltp = self.dict['minAltP']                                                           
		
                # Reading in the DEM
		if self.verb:
			print 'PROGRESS: READING DEM'
                fin = open(self.hfile,'rb');
                if self.fmt in ('HGT'):
                    dem = np.fromfile(file=fin,dtype=self.demtype,count=self.nx*self.ny).reshape(self.ny,self.nx)
                elif self.fmt in ('RMG'):
                    dem = np.fromfile(file=fin,dtype=self.demtype,count=2*self.nx*self.ny).reshape(self.ny,2*self.nx)
                    dem = dem[:,self.nx:]
                dem = np.round(dem).astype(np.int)
		fin.close()

		# check output, and open file if necessary
		outFile = isinstance(dataobj,str)
                if outFile:
                        fout = open(dataobj,'wb')
                        dout = np.zeros((self.ny,self.nx))
                else:
                        assert ((dataobj.shape[0]==self.ny) & (dataobj.shape[1]==self.nx)), 'PyAPS: Not a valid data object.'
                        dout = dataobj

		# Incidence
                cinc = np.cos(inc*np.pi/180.0)

		# Create the 1d interpolator
		if self.verb:
			print 'PROGRESS: FINE INTERPOLATION OF HEIGHT LEVELS'
                intp_1d = si.interp1d(self.hgt,self.Delfn,kind='cubic',axis=1)

		# Interpolate the delay function every meter
		dem[dem < minAltp] = minAltp
                minH = dem.min()
                maxH = dem.max()+1
		kh = np.arange(minH,maxH)
		Delfn_1m = intp_1d(kh)
		self.Delfn_interp = Delfn_1m.copy()
		self.alti = kh	

		# Reshape Delfn
		Lonu = np.unique(self.lonlist)
		Latu = np.unique(self.latlist)
		nLon = len(Lonu)
		nLat = len(Latu)
		Delfn_1m = np.reshape(Delfn_1m.T,(len(kh),nLat,nLon))
		self.Delfn_1m = Delfn_1m
		
		# build the x array
		xarr = np.arange(1.,self.nx+1.)

		# Create the cube interpolator for the bilinear method
		if self.verb:
			print 'PROGRESS: CREATE THE BILINEAR INTERPOLATION FUNCTION'
		bilicube = processor.Bilinear2DInterpolator(Lonu,Latu,Delfn_1m,cube=True)

		# Loop on the lines
		if self.verb:
			toto = utils.ProgressBar(maxValue=self.ny)
			print 'PROGRESS: MAPPING THE DELAY'
		for m in xrange(self.ny):
			if self.verb:
				toto.update(m,every=5)

			# Transfert (range,azimuth) to (lon,lat)
			yarr = (m+1)*np.ones((xarr.shape))
			[loni,lati] = utils.rdr2glob(self.nx,self.ny,self.lat,self.lon,xarr,yarr)

			# Make the bilinear interpolation
			D = dem[m,:] - minH
			val = bilicube(loni,lati,D)*np.pi*4.0/(cinc*wvl)
			
			if outFile:
				resy = val.astype(np.float32)
				resy.tofile(fout)
			else:
				dataobj[m,:] = val

		if self.verb:
			toto.close()

		# Close if outfile		
		if outFile:
                        fout.close()

	def getgeodelay(self,dataobj,lat=None,lon=None,inc=None,wvl=np.pi*4.):
                '''Write delay to a matrix / HDF5 object or a file directly. This is used when the latitude and longitude values are available for each radar pixel. Incidence angle can be a constant or a file name.
		Bilinear Interpolation is used.
                
                Args:   
                        * dataobj  (str or HDF5 or np.array): Final output. If str, output is written to file.
                Kwargs:         
                        * lat  (str)            : Path to the latitude file (np.float32).
                        * lon  (str)            : Path to the longitude file (np.float32).
                        * inc  (str or np.float): Path to incidence angle file in degrees (str) or a constant float.
                        * wvl  (np.float)       : Wavelength in meters. Default output results in delay in meters.
                
                .. note::
                        If dataobj is string, output is written to the file.
                        If np.array or HDF5 object, it should be of size (ny,nx).'''
	
		assert lat is not None, 'PyAPS: Need a valid latitude file.'
                assert lon is not None, 'PyAPS: Need a valid longitude file.'

                if isinstance(inc,float) or isinstance(inc,np.float64) or isinstance(inc,np.float32):
                        cinc = np.cos(inc*np.pi/180.)
                        incFileFlag = False
                else:
                        assert inc is not None, 'PyAPS: Need a valid incidence angle file or constant.'
                        incFileFlag = True
                        incin = open(inc,'rb')
                
                latin = open(lat,'rb')
                lonin = open(lon,'rb')

		minAltp = self.dict['minAltP']

                # Reading in the DEM
		if self.verb:
	                print 'PROGRESS: READING DEM'
                fin = open(self.hfile,'rb');
                if self.fmt in ('HGT'):
                    dem = np.fromfile(file=fin,dtype=self.demtype,count=self.nx*self.ny).reshape(self.ny,self.nx)
                elif self.fmt in ('RMG'):
                    dem = np.fromfile(file=fin,dtype=self.demtype,count=2*self.nx*self.ny).reshape(self.ny,2*self.nx)
                    dem = dem[:,self.nx:]
                dem = np.round(dem).astype(np.int)
                fin.close()

                # check output, and open file if necessary
                outFile = isinstance(dataobj,str)
                if outFile:
                        fout = open(dataobj,'wb')
                        dout = np.zeros((self.ny,self.nx))
                else:
                        assert ((dataobj.shape[0]==self.ny) & (dataobj.shape[1]==self.nx)), 'PyAPS: Not a valid data object.'
                        dout = dataobj

		# Create the 1d interpolator
		if self.verb:
        	        print 'PROGRESS: FINE INTERPOLATION OF HEIGHT LEVELS'
                intp_1d = si.interp1d(self.hgt,self.Delfn,kind='cubic',axis=1)

                # Interpolate the delay function every meter
                dem[dem < minAltp] = minAltp
                minH = dem.min()
                maxH = dem.max()+1
                kh = np.arange(minH,maxH)
                Delfn_1m = intp_1d(kh)
                self.Delfn_interp = Delfn_1m.copy()
                self.alti = kh

		# Reshape Delfn
                Lonu = np.unique(self.lonlist)
                Latu = np.unique(self.latlist)
                nLon = len(Lonu)
                nLat = len(Latu)
                Delfn_1m = np.reshape(Delfn_1m.T,(len(kh),nLat,nLon))
                self.Delfn_1m = Delfn_1m

                # build the x array
                xarr = np.arange(1.,self.nx+1.)

                # Create the cube interpolator for the bilinear method
		if self.verb:
	                print 'PROGRESS: CREATE THE BILINEAR INTERPOLATION FUNCTION'
                bilicube = processor.Bilinear2DInterpolator(Lonu,Latu,Delfn_1m,cube=True)

		# Loop on the lines
		if self.verb:
	                toto = utils.ProgressBar(maxValue=self.ny)
       		        print 'PROGRESS: MAPPING THE DELAY'
                for m in xrange(self.ny):
			if self.verb:
	                        toto.update(m,every=5)
		
			# Get Latitude and Longitude arrays
			lati = np.fromfile(file=latin,dtype=np.float32,count=self.nx)
                        loni = np.fromfile(file=lonin,dtype=np.float32,count=self.nx)
			loni[loni<0.] += 360.
			ii = np.where(np.isnan(lati))
			jj = np.where(np.isnan(loni))
			xx = np.union1d(ii,jj)
			lati[xx]=0.0
			loni[xx]=0.0

			# Get incidence if file provided
			if incFileFlag:
				incz = np.fromfile(file=incin,dtype=np.float32,count=self.nx)
				cinc = np.cos(incz*np.pi/180.)

			# Make the Bilinear interpolation
			D = dem[m,:] - minH
                        val = bilicube(loni,lati,D)*np.pi*4.0/(cinc*wvl)
			val[xx] = np.nan

                        if outFile:
                                resy = val.astype(np.float32)
                                resy.tofile(fout)
                        else:
                                dataobj[m,:] = val

		if self.verb:
			toto.close()
		latin.close()
		lonin.close()
		if incFileFlag:
			incin.close()

		if outFile:
			fout.close()

	def geomerisfactor(self,dataobj,lat=None,lon=None,inc=None,wvl=np.pi*4.):
		''' Write pi-factor from Li et al 2012 to a matrix / HDF5 object or a file directly. This is used when the latitude and longitude values are available for each radar pixel. Incidence angle can be a constant or a file name.

		Bilinear Interpolation is used.

		Args:
			* dataobj (str or HDF5 or np.array): Final output. If str, output is written to file.
                Kwargs:         
                        * lat  (str)            : Path to the latitude file (np.float32).
                        * lon  (str)            : Path to the longitude file (np.float32).
                        * inc  (str or np.float): Path to incidence angle file in degrees (str) or a constant float.
                        * wvl  (np.float)       : Wavelength in meters. Default output results in delay in meters.
                
                .. note::
                        If dataobj is string, output is written to the file.
                        If np.array or HDF5 object, it should be of size (ny,nx).'''

		assert lat is not None, 'PyAPS: Need a valid latitude file.'
                assert lon is not None, 'PyAPS: Need a valid longitude file.'

		if isinstance(inc,float) or isinstance(inc,np.float64) or isinstance(inc,np.float32):
                        cinc = np.cos(inc*np.pi/180.)
                        incFileFlag = False
                else:
                        assert inc is not None, 'PyAPS: Need a valid incidence angle file or constant.'
                        incFileFlag = True
                        incin = open(inc,'rb')

                latin = open(lat,'rb')
                lonin = open(lon,'rb')

                minAltp = self.dict['minAltP']

		# Compute the two integrals
		WonT = self.Vi/self.Ti
		WonT2 = WonT/self.Ti

		S1 = intg.cumtrapz(WonT,x=self.hgt,axis=-1)
		val = 2*S1[:,-1]-S1[:,-2]
	        val = val[:,None]
	        S1 = np.concatenate((S1,val),axis=-1)
		del WonT
                S2 = intg.cumtrapz(WonT2,x=self.hgt,axis=-1)
                val = 2*S2[:,-1]-S2[:,-2]
                val = val[:,None]
                S2 = np.concatenate((S2,val),axis=-1)
                del WonT2
		Tm = S1/S2
		self.Tm = Tm
	
                # Reading in the DEM
                if self.verb:
                        print 'PROGRESS: READING DEM'
                fin = open(self.hfile,'rb');
                if self.fmt in ('HGT'):
                    dem = np.fromfile(file=fin,dtype=self.demtype,count=self.nx*self.ny).reshape(self.ny,self.nx)
                elif self.fmt in ('RMG'):
                    dem = np.fromfile(file=fin,dtype=self.demtype,count=2*self.nx*self.ny).reshape(self.ny,2*self.nx)
                    dem = dem[:,self.nx:]
                dem = np.round(dem).astype(np.int)
                fin.close()
		
                # check output, and open file if necessary
                outFile = isinstance(dataobj,str)
                if outFile:
                        fout = open(dataobj,'wb')
                        dout = np.zeros((self.ny,self.nx))
                else:
                        assert ((dataobj.shape[0]==self.ny) & (dataobj.shape[1]==self.nx)), 'PyAPS: Not a valid data object.'
                        dout = dataobj

		# Create the 1d interpolator
                if self.verb:
                        print 'PROGRESS: FINE INTERPOLATION OF HEIGHT LEVELS'
                intp_1d = si.interp1d(self.hgt,Tm,kind='cubic',axis=1)

		# Interpolate the Tm variable every meter
                dem[dem < minAltp] = minAltp
                minH = dem.min()
                maxH = dem.max()+1
                kh = np.arange(minH,maxH)
                Tm_1m = intp_1d(kh)
                self.alti = kh

		# Reshape Tm
                Lonu = np.unique(self.lonlist)
                Latu = np.unique(self.latlist)
                nLon = len(Lonu)
                nLat = len(Latu)
                Tm_1m = np.reshape(Tm_1m.T,(len(kh),nLat,nLon))
		self.Tm_1m = Tm_1m

		# build the x array
                xarr = np.arange(1.,self.nx+1.)

		# Create the cube interpolator for the bilinear method
                if self.verb:   
                        print 'PROGRESS: CREATE THE BILINEAR INTERPOLATION FUNCTION'
                bilicube = processor.Bilinear2DInterpolator(Lonu,Latu,Tm_1m,cube=True)

		# Get the values from the dictionnary
		k1 = self.dict['k1']
		k2 = self.dict['k2']
                k3 = self.dict['k3']
	        mmO = self.dict['mmO']
	        mmH = self.dict['mmH']
		mma = self.dict['mma']
		w = (2*mmH + mmO)/mma
		Rv = self.dict['Rv']
		Rho = self.dict['Rho']

                # Loop on the lines
                if self.verb:
                        toto = utils.ProgressBar(maxValue=self.ny)
                        print 'PROGRESS: MAPPING THE DELAY'
                for m in xrange(self.ny):
                        if self.verb:
                                toto.update(m,every=5)
			# Get Latitude and Longitude arrays
                        lati = np.fromfile(file=latin,dtype=np.float32,count=self.nx)
                        loni = np.fromfile(file=lonin,dtype=np.float32,count=self.nx)
                        loni[loni<0.] += 360.
                        ii = np.where(np.isnan(lati))
                        jj = np.where(np.isnan(loni))
                        xx = np.union1d(ii,jj)
                        lati[xx]=0.0
                        loni[xx]=0.0

			# Get incidence if file provided
                        if incFileFlag:
                                incz = np.fromfile(file=incin,dtype=np.float32,count=self.nx)
                                cinc = np.cos(incz*np.pi/180.)
		
			# Make the Bilinear interpolation
                        D = dem[m,:] - minH
                        val = bilicube(loni,lati,D)
			val = 0.000001 * Rho * Rv * ( k3/val + k2 - w*k1) * np.pi*4.0/(cinc*wvl)
                        val[xx] = np.nan

			if outFile:
                                resy = val.astype(np.float32)
                                resy.tofile(fout)
                        else:
                                dataobj[m,:] = val

                if self.verb:
                        toto.close()
                latin.close()
                lonin.close()
                if incFileFlag:
                        incin.close()

                if outFile:
                        fout.close()

	def merisfactor(self,dataobj,inc=0.0,wvl=np.pi*4.):
		''' Write pi-factor from Li et al 2012 to a matrix / HDF5 object or a file directly. This is used when no lat/lon file is known.

                Bilinear Interpolation is used.

                Args:
                        * dataobj (str or HDF5 or np.array): Final output. If str, output is written to file.
                Kwargs:         
                        * inc  (np.float)       : Incidence angle in degrees.
                        * wvl  (np.float)       : Wavelength in meters. Default output results in delay in meters.
                
                .. note::
                        If dataobj is string, output is written to the file.
                        If np.array or HDF5 object, it should be of size (ny,nx).'''

		
		minAltp = self.dict['minAltP']

                # Incidence
                cinc = np.cos(inc*np.pi/180.0)

		# Compute the two integrals
		WonT = self.Vi/self.Ti
		WonT2 = WonT/self.Ti

		S1 = intg.cumtrapz(WonT,x=self.hgt,axis=-1)
		val = 2*S1[:,-1]-S1[:,-2]	
		val = val[:,None]
		S1 = np.concatenate((S1,val),axis=-1)
		del WonT   
		S2 = intg.cumtrapz(WonT2,x=self.hgt,axis=-1)
		val = 2*S2[:,-1]-S2[:,-2] 	
		val = val[:,None] 
		S2 = np.concatenate((S2,val),axis=-1) 
		del WonT2 
		Tm = S1/S2 
		self.Tm = Tm                                                                                                        
	
                # Reading in the DEM                                                                                                 
                if self.verb:                                                                                                        
                        print 'PROGRESS: READING DEM'                                                                                
                fin = open(self.hfile,'rb');                                                                                         
                if self.fmt in ('HGT'):
                    dem = np.fromfile(file=fin,dtype=self.demtype,count=self.nx*self.ny).reshape(self.ny,self.nx)                    
                elif self.fmt in ('RMG'):
                    dem = np.fromfile(file=fin,dtype=self.demtype,count=2*self.nx*self.ny).reshape(self.ny,2*self.nx)                
                    dem = dem[:,self.nx:]                                                                                            
                dem = np.round(dem).astype(np.int)
                fin.close()
                
                # check output, and open file if necessary                                                                           
                outFile = isinstance(dataobj,str)                                                                                    
                if outFile:
                        fout = open(dataobj,'wb')
                        dout = np.zeros((self.ny,self.nx))                                                                           
                else:
                        assert ((dataobj.shape[0]==self.ny) & (dataobj.shape[1]==self.nx)), 'PyAPS: Not a valid data object.'
                        dout = dataobj

                # Create the 1d interpolator
                if self.verb:
                        print 'PROGRESS: FINE INTERPOLATION OF HEIGHT LEVELS'
                intp_1d = si.interp1d(self.hgt,Tm,kind='cubic',axis=1)

                # Interpolate the Tm variable every meter
                dem[dem < minAltp] = minAltp
                minH = dem.min()
                maxH = dem.max()+1
                kh = np.arange(minH,maxH)
                Tm_1m = intp_1d(kh)
                self.alti = kh

                # Reshape Tm
                Lonu = np.unique(self.lonlist)
                Latu = np.unique(self.latlist)
                nLon = len(Lonu)
                nLat = len(Latu)
                Tm_1m = np.reshape(Tm_1m.T,(len(kh),nLat,nLon))
                self.Tm_1m = Tm_1m

                # build the x array
                xarr = np.arange(1.,self.nx+1.)

		# Create the cube interpolator for the bilinear method
                if self.verb:
                        print 'PROGRESS: CREATE THE BILINEAR INTERPOLATION FUNCTION'
                bilicube = processor.Bilinear2DInterpolator(Lonu,Latu,Tm_1m,cube=True)

		# Get the values from the dictionnary
                k1 = self.dict['k1']
                k2 = self.dict['k2']
                k3 = self.dict['k3']
                mmO = self.dict['mmO']
                mmH = self.dict['mmH']
                mma = self.dict['mma']
                w = (2*mmH + mmO)/mma
                Rv = self.dict['Rv']
                Rho = self.dict['Rho']

                # Loop on the lines
                if self.verb:
                        toto = utils.ProgressBar(maxValue=self.ny)
                        print 'PROGRESS: MAPPING THE DELAY'
                for m in xrange(self.ny):
                        if self.verb:
                                toto.update(m,every=5)

                        # Transfert (range,azimuth) to (lon,lat)
                        yarr = (m+1)*np.ones((xarr.shape))
                        [loni,lati] = utils.rdr2glob(self.nx,self.ny,self.lat,self.lon,xarr,yarr)

                        # Make the bilinear interpolation
                        D = dem[m,:] - minH
                        val = bilicube(loni,lati,D)
			val = 0.000001 * Rho * Rv * ( k3/val + k2 - w*k1) * np.pi*4.0/(cinc*wvl)

                        if outFile:
                                resy = val.astype(np.float32)
                                resy.tofile(fout)
                        else:
                                dataobj[m,:] = val

                if self.verb:
                        toto.close()

                # Close if outfile              
                if outFile:
                        fout.close()

###########End of PyAPS_rdr class###############################



############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
############################################################
