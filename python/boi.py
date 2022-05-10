from scipy.constants import speed_of_light
from scipy.constants import pi

def s1_azfm(r, t0, azp):
  """azfr = s1_azfm(r, t0, azp)
  Calculate azimuth FM rate given slant range, reference slant-range delay and the azimuth FM rate polynomial for ScanSAR data
  **Arguments:**
  * r:    slant range (meters)
  * t0:   reference slant range time for the polynomial (center swath delay in s)
  * azp:  polynomial coefficients
  **Output:**
  * the function returns the azimuth FM rate"""
  tsr = 2.0 * r / speed_of_light
  dt = tsr - t0
  azfr = azp[0] + dt * (azp[1] + dt*(azp[2] + dt*(azp[3] + dt*azp[4])))
  return azfr


def get_param_gamma(param, parfile, floatt = True, pos = 0):
    a = grep1line(param,parfile).split()[1+pos]
    if floatt:
        a = float(a)
    return a


def get_dfDC(path_to_slcdir, f0=5405000500, burst_interval = 2.758277, returnka = True, returnperswath = False):
    #f0 = get_param_gamma('radar_frequency', parfile)
    #burst_interval = get_param_gamma('burst_interval', topsparfile)
    parfile = glob.glob(path_to_slcdir+'/????????.slc.par')[0]
    topsparfiles = glob.glob(path_to_slcdir+'/????????.IW?.slc.TOPS_par')
    iwparfiles = glob.glob(path_to_slcdir+'/????????.IW?.slc.par')
    #
    lam = speed_of_light / f0
    dfDC = []
    kas = []
    #krs = []
    #print('This is a proper solution but applied to primary SLC image. originally it is applied by GAMMA on the RSLC...')
    for n in range(len(topsparfiles)):
        topsparfile = topsparfiles[n]
        iwparfile = iwparfiles[n]
        az_steering_rate = get_param_gamma('az_steering_rate', topsparfile) # az_steering_rate is the antenna beam steering rate
        r1 = get_param_gamma('center_range_slc', iwparfile)
        #get the satellite velocity
        #midNstate = int(get_param_gamma('number_of_state_vectors', iwparfile)/2)+1
        # ... actually number of burst info differs... so just using the 1st burst - as anyway we do quite drastic change to dfDC - mean from swaths
        midNstate = 1
        sv = 'state_vector_velocity_' + str(midNstate)
        velc1 = get_param_gamma(sv, iwparfile, pos=0)
        velc2 = get_param_gamma(sv, iwparfile, pos=1)
        velc3 = get_param_gamma(sv, iwparfile, pos=2)
        vsat = np.sqrt(velc1**2 + velc2**2 + velc3**2)
        # now some calculations
        afmrate_srdelay = get_param_gamma('az_fmrate_srdelay_'+ str(midNstate), topsparfile)
        afmrate_poly = []
        afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 0))
        afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 1))
        afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 2))
        try:
            afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 3))
        except:
            afmrate_poly.append(0)
        try:
            afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 4))
        except:
            afmrate_poly.append(0)
        ka = s1_azfm(r1, afmrate_srdelay, afmrate_poly) #unit: Hz/s == 1/s^2
        kr = -2.0 * vsat * az_steering_rate*(pi / 180.0) / lam
        if (kr != 0.0):
            #kt = ka * kr/(kr - ka)
            # but maybe should be kt = (kr*ka)/(ka-kr) # see https://iopscience.iop.org/article/10.1088/1755-1315/57/1/012019/pdf  --- and remotesensing-12-01189-v2, and Fattahi et al...
            # ok, gamma reads kr to be ka... corrected
            kt = kr * ka/(ka - kr)
        else:
            kt = -ka
        #finally calculate dfDC:
        #burst_interval = get_param_gamma('burst_interval', topsparfile)
        kas.append(ka)
        #krs.append(kr)
        dfDC.append(kt*burst_interval) #burst_interval is time within the burst... we can also just calculate.. see Grandin: eq 15: hal.archives-ouvertes.fr/hal-01621519/document
        #ok, that's the thing - burst_interval is actually t(n+1) - t(n) - see remotesensing-12-01189-v2
        #so it should be kt * -burst_interval, that is why GAMMA has the -kt J ... ok, good to realise this
    if not returnperswath:
        dfDC = np.mean(dfDC)
        ka = np.mean(kas)
    #kr = np.mean(krs)
    if returnka:
        return dfDC, ka #, kr
    else:
        return dfDC


'''
# extract overlaps and make their ifg and coh
ScanSAR_burst_overlap tab/20220309R_tab 20220309_ovl 4 20 - - tab/20190301R_tab
ScanSAR_burst_overlap tab/20220321R_tab 20220321_ovl 4 20 - - tab/20190301R_tab
create_offset 20220309_ovl.bwd.slc.par 20220321_ovl.bwd.slc.par ifg.off 1 20 4 0
phase_sim_orb 20220309_ovl.bwd.slc.par 20220321_ovl.bwd.slc.par ifg.off geo/20190301.hgt 20220309_20220321.simunw RSLC/20190301/20190301.rslc.par 
SLC_diff_intf 20220309_ovl.bwd.slc 20220321_ovl.bwd.slc 20220309_ovl.bwd.slc.par 20220321_ovl.bwd.slc.par ifg.off 20220309_20220321.simunw bwd.ifg 20 4 0 0
SLC_diff_intf 20220309_ovl.fwd.slc 20220321_ovl.fwd.slc 20220309_ovl.fwd.slc.par 20220321_ovl.fwd.slc.par ifg.off 20220309_20220321.simunw fwd.ifg 20 4 0 0
python3 -c "import numpy as np; a=np.fromfile('fwd.ifg', np.complex64).byteswap(); b=np.fromfile('bwd.ifg', np.complex64).byteswap(); c = a*np.conj(b);  c.byteswap().tofile('ddiff')"
# hopefully ok for coh...:
cc_wave ddiff RSLC/20220309/20220309.rslc.mli RSLC/20220321/20220321.rslc.mli ddiff_coh 3388
cpx_to_real ddiff ddiff.pha 3388 4
geocode_back ddiff_coh 3388 geo/20190301.lt_fine ddiff_coh.geo 2548 1250 0 0
geocode_back ddiff.pha 3388 geo/20190301.lt_fine ddiff.pha.geo 2548 1250 0 0
data2geotiff geo/EQA.dem_par ddiff_coh.geo 2 ddiff_coh.geo.tif 0.0
data2geotiff geo/EQA.dem_par ddiff.pha.geo 2 ddiff.geo.tif 0.0

### not this below anymore
cpx_to_real bwd.ifg bwd.ifg.pha 3388 4
cpx_to_real fwd.ifg fwd.ifg.pha 3388 4
#geocode_back bwd.ifg 3388 geo/20190301.lt_fine bwd.ifg.geo 2548 1250 1 1
#geocode_back fwd.ifg 3388 geo/20190301.lt_fine fwd.ifg.geo 2548 1250 1 1
geocode_back bwd.ifg.pha 3388 geo/20190301.lt_fine bwd.ifg_pha.geo 2548 1250 0 0
geocode_back fwd.ifg.pha 3388 geo/20190301.lt_fine fwd.ifg_pha.geo 2548 1250 0 0
data2geotiff geo/EQA.dem_par bwd.ifg_pha.geo 2 bwd.ifg_pha.geo.tif 0.0
data2geotiff geo/EQA.dem_par fwd.ifg_pha.geo 2 fwd.ifg_pha.geo.tif 0.0
gmt grdimage fwd.ifg_pha.geo.tif -C$LiCSARpath/misc/pha.cpt -JM1 -nn+t0.1 -Q -Afwd.pha.png
gmt grdimage bwd.ifg_pha.geo.tif -C$LiCSARpath/misc/pha.cpt -JM1 -nn+t0.1 -Q -Abwd.pha.png

'''

# load to xr and do difference, and recomputation to N-S displacements there:
from unwrp_multiscale import *
incoh = load_tif2xr('ddiff_coh.geo.tif')
inpha = load_tif2xr('ddiff.geo.tif')
ifg = xr.Dataset()
ifg['pha'] = inpha
ifg['coh'] = ifg['pha']
ifg['coh'].values = incoh.val
inmask = incoh.copy(deep=True)
inmask.values = np.byte(incoh > 0)
ifg['mask'] = ifg['pha']
ifg['mask'].values = inmask.values
ifg['cpx'] = ifg.coh.copy()
ifg['cpx'].values = magpha2RI_array(ifg.coh.values, ifg.pha.values)

ref = ifg.sel(lon=slice(-28,-27.5), lat=slice(38.58, 38.5))
# get overall multilooked phase, as reference:
refavg = ref.cpx.sum()/ref.cpx.count()
refavg = np.angle(refavg)

aoi = ifg.sel(lon=slice(-28.4,-28.2), lat=slice(38.504, 38.444))
aoi['pha'].values = wrap2phase(aoi.pha - refavg)
aoi = filter_ifg_ml(aoi)
aoi['cpx'].values = magpha2RI_array(aoi.coh.values, aoi.gauss_pha.values)

aoiavg = aoi.cpx.sum()/aoi.cpx.count()
aoiavg = np.angle(aoiavg)

#using average dfDC and resolution, max diff between swaths 1 and 3 is approx 1rad = 14mm, so mean would differ max 1 rad = 7 mm
path_to_slcdir = 'SLC/20190301'
dfdc = get_dfDC(path_to_slcdir,returnka = False)
parfile = glob.glob(path_to_slcdir+'/????????.slc.par')[0]
azires = get_param_gamma('azimuth_pixel_spacing',parfile)
prf = 486.486
rad2mm_coef = -prf*azires*1000 / (2*pi*dfdc)







ddiff = 

fwd = load_tif2xr('fwd.ifg_pha.geo.tif')
bwd = load_tif2xr('bwd.ifg_pha.geo.tif')
diff=fwd.copy(deep=True).rename('diff')
diff.values = wrap2phase(fwd - bwd)
diff = diff.where(diff != 0)
# refer to ref point
reflon=-27.76415
reflat=38.54613
diff = diff - diff.sel(lon=reflon, lat=reflat, method='nearest')
diff.values = wrap2phase(diff.values)
# so now do some more to unwrap:
#mag = np.ones(diff.shape)
coh = 
cpx = magpha2RI_array(coh.values,diff.values)
# refarea:
refdiff = diff.sel(lon=slice(-28,-27.5), lat=slice(38.58, 38.5))
refcpx = refdiff.copy(deep=True)
mag = np.ones(refdiff.shape)
refcpx.values = magpha2RI_array(mag,refdiff.values)
# get overall multilooked phase, as reference:
avg = refcpx.sum()/refcpx.count()
avg = np.angle(avg)
# aoi
#refdiff.values = wrap2phase(refdiff.values - avg)

aoidiff = diff.sel(lon=slice(-28.4,-28.2), lat=slice(38.504, 38.444))
aoidiff.values = wrap2phase(aoidiff - avg)
unw = 

#using average dfDC and resolution, max diff between swaths 1 and 3 is approx 1rad = 14mm, so mean would differ max 1 rad = 7 mm
path_to_slcdir = 'SLC/20190301'
dfdc = get_dfDC(path_to_slcdir,returnka = False)
parfile = glob.glob(path_to_slcdir+'/????????.slc.par')[0]
azires = get_param_gamma('azimuth_pixel_spacing',parfile)
prf = 486.486
rad2mm_coef = -prf*azires*1000 / (2*pi*dfdc)



diffmm = diff*rad2mm_coef #*14000 
38.44463577, -28.40075818
38.50442392, -28.20322605
aoidiffmm = aoidiff*rad2mm_coef


refcpx = cpx.where(lon<-28).where(lat<38.57).where(
cpx.sel(lon=slice(-28,-27.5), lat=slice(38.42, 38.58))
# here maybe unwrap diff?
diffpx = diff*coef
diffmm = diffpx*14000 # approx.
diffmm.rename('mm').plot(); plt.show()
