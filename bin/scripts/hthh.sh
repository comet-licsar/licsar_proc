# this is to process HTHH frame - useful later to do crop etc. in LiCSAR...

cd $BATCH_CACHE_DIR/102D_11050_000002

# 1. init
if [ 1 == 0 ]; then
#to create a crop
SLC_mosaic_S1_TOPS 20211210_tab 20211210.slc 20211210.slc.par 1 1 
SLC_copy 20211210.slc 20211210.slc.par 20211210.crop.slc 20211210.crop.slc.par - - 4500 1300 1120 500
multi_look 20211210.crop.slc 20211210.crop.slc.par 20211210.crop.mli 20211210.crop.mli.par 1 1
mliwid=1300
#cpxfiddle -w $mliwid -o sunraster -c gray -q normal -f r4 -M1/1 -Bb 20211210.crop.mli > 20211210.crop.mli.ras

#to DEM-geocode:
gc_map2 $mlipar geo_crop/dem_crop.dem_par DEM/dem_crop.dem geo_crop/dem_seg_par geo_crop/dem_seg geo_crop/lt 0.5 0.5 geo_crop/lsmap geo_crop/lsmap_rdc geo_crop/inc geo_crop/res geo_crop/offnadir geo_crop/sim_sar geo_crop/u geo_crop/v geo_crop/psi geo_crop/pix - - - 12 - - 

#check output!
create_diff_par $mlipar - geo_crop/$m.diff_par 1
init_offsetm $mli geo/sim_sar geo_crop/$m.diff_par

#to gec-geocode:
cd GEC
create_dem_par fake.dem_par $mlipar - - - 4326 0
gec_map $slc.par - fake.dem_par 0 fake.seg.dem_par gec.lt - -
#width: 226   lines: 275
demwid=226
geocode_back $mli $mliwid gec.lt mli.geo $demwid  
dispwr mli.geo $demwid
fi


# 2. crop all slcs, and compute offsets, separately for each island:
m=20211210
slcdir=/work/scratch-pw/earmla/LiCSAR_temp/batchdir/102D_11050_000002/SLC
rslcdir=/work/scratch-pw/earmla/LiCSAR_temp/batchdir/102D_11050_000002/RSLC

if [ 1 == 1 ]; then
# crop all SLCs
cd SLC
for iscode in full L R; do

if [ $iscode == 'full' ]; then
#for whole island area:
#iscode='full'
sr=4500
nr=1300
sa=1120
na=500
elif [ $iscode == 'L' ]; then
#for left island, i.e. the one further in range:
#iscode='L'
let sr=4500+889
let nr=254
let sa=1120+209
let na=174
elif [ $iscode == 'R' ]; then
#for right island:
#iscode='R'
let sr=4500+140
let nr=558
let sa=1120+135
let na=72
fi

rgofffile=/work/scratch-pw/earmla/LiCSAR_temp/batchdir/102D_11050_000002/rgoffsets.$iscode.csv
echo "epoch,rg_orb,rg_off" > $rgofffile
echo $m",0,0" >> $rgofffile

for x in `ls`; do
 cd $x
 echo "doing "$x
 if [ ! -f $x.der.slc ]; then
  echo "deramping and mosaicking"
  echo $x.IW3.slc $x.IW3.slc.par $x.IW3.slc.TOPS_par > $x.tab
  echo $x.IW3.der.slc $x.IW3.der.slc.par $x.IW3.der.slc.TOPS_par > $x.der.tab
  SLC_deramp_S1_TOPS $x.tab $x.der.tab 0 1
  #SLC_mosaic_S1_TOPS $x.tab $x.slc $x.slc.par 1 1 0
  SLC_mosaic_S1_TOPS $x.der.tab $x.der.slc $x.der.slc.par 1 1 0
 fi
 if [ ! -f $x.der.crop.$iscode.mli.png ]; then
     echo "cropping"
     #SLC_copy $x.slc $x.slc.par $x.crop.slc $x.crop.slc.par - - 4500 $nr 1120 500
     SLC_copy $x.der.slc $x.der.slc.par $x.der.crop.$iscode.slc $x.der.crop.$iscode.slc.par - - $sr $nr $sa $na
     #multi_look $x.crop.slc $x.crop.slc.par $x.crop.mli $x.crop.mli.par 1 1
     multi_look $x.der.crop.$iscode.slc $x.der.crop.$iscode.slc.par $x.der.crop.$iscode.mli $x.der.crop.$iscode.mli.par 1 1
     raspwr $x.der.crop.$iscode.mli $nr
     convert $x.der.crop.$iscode.mli.bmp $x.der.crop.$iscode.mli.png
     #display $x.der.crop.$iscode.mli.png
 fi
 cd ..
done

for s in `ls`; do
 cd $s
 echo "doing "$s
 # now estimating offsets
 if [ ! $m == $s ]; then
  ext=der.crop.$iscode
  mpar=$slcdir/$m/$m.$ext.slc.par
  mslc=$slcdir/$m/$m.$ext.slc
  spar=$slcdir/$s/$s.$ext.slc.par
  sslc=$slcdir/$s/$s.$ext.slc
  mkdir -p $rslcdir/$s
  rslc=$rslcdir/$s/$s.$ext.rslc
  off=$rslcdir/$s/$m'_'$s.$ext.off
  # use orbits only
  create_offset $mpar $spar $off 1 1 1 0
  init_offset_orbit $mpar $spar $off
  cp $off $off.orbonly
  # now do CC
  #offset_pwr $mslc $sslc $mpar $spar $off $off.offs $off.ccp 64 64 $off.offsets 1 16 16 0.1
  #now fit the offsets and improve in off file, if possible.. use only polynomial 1 (translation)!
  #offset_fit $m'_'$s.offs ccp $m'_'$s.off coffs coffsets 0.15 1
  # refinement
  offset_pwr $mslc $sslc $mpar $spar $off $off.offs $off.ccp 32 32 $off.offsets 1 32 32 0.5
  offset_fit $off.offs $off.ccp $off $off.coffs $off.coffsets 0.5 1
  offset_pwr $mslc $sslc $mpar $spar $off $off.offs $off.ccp 16 16 $off.offsets 2 32 16 0.5
  offset_fit $off.offs $off.ccp $off $off.coffs $off.coffsets 0.6 1
  #offset_pwr $mslc $sslc $mpar $spar $off $off.offs $off.ccp 16 16 $off.offsets 2 32 32 0.5
  # now extract the values:
  rgorb=`grep range_offset_polynomial $off.orbonly | gawk {'print $2'}`
  rgoff=`grep range_offset_polynomial $off | gawk {'print $2'}`
  echo $s','$rgorb','$rgoff >> $rgofffile
  if [ $iscode == 'full' ]; then
   # 3b. coregister the mosaics, but only using the full area..
   # interpolate to RSLC
   SLC_interp $sslc $mpar $spar $off $rslc $rslc.par - -
  fi
 fi
 cd ..
done

cd ..
fi

done



# 2. get tides:
out_SET=tides.txt
echo "epoch,E,N,U" > $out_SET
lon=-175.3913
lat=-20.54436
masterdate=20211210
masterdate=`echo ${masterdate:0:4}-${masterdate:4:2}-${masterdate:6:2}`
centertime="17:09:01"
masterdt=$masterdate"T"$centertime
#heading=193
mtide=`gmt earthtide -L$lon/$lat -T$masterdt 2>/dev/null | sed 's/\t/,/g'`
NM=`echo $mtide | cut -d ',' -f2`
EM=`echo $mtide | cut -d ',' -f3`
VM=`echo $mtide | cut -d ',' -f4`

for epochdate in `ls RSLC`; do
 echo $epochdate
 epochdatee=`echo ${epochdate:0:4}-${epochdate:4:2}-${epochdate:6:2}`
 epochdt=$epochdatee"T"$centertime
 etide=`gmt earthtide -L$lon/$lat -T$epochdt 2>/dev/null | sed 's/\t/,/g'`
 NE=`echo $etide | cut -d ',' -f2`
 EE=`echo $etide | cut -d ',' -f3`
 VE=`echo $etide | cut -d ',' -f4`
 U=`python3 -c "print(("$VE")-("$VM"))"`
 E=`python3 -c "print(("$EE")-("$EM"))"`
 N=`python3 -c "print(("$NE")-("$NM"))"`
 echo $epochdate","$E","$N","$U >> $out_SET
done



# 3. now python - load it and plot it babe!
import matplotlib.pyplot as plt

import pandas as pd
L=pd.read_csv('rgoffsets.L.csv')
R=pd.read_csv('rgoffsets.R.csv')
full=pd.read_csv('rgoffsets.full.csv')
tides=pd.read_csv('tides.txt')
L = L.set_index(pd.to_datetime(L['epoch'].astype(str)))
R = R.set_index(pd.to_datetime(R['epoch'].astype(str)))
tides = tides.set_index(pd.to_datetime(tides['epoch'].astype(str)))
full = full.set_index(pd.to_datetime(full['epoch'].astype(str)))
R = R.query('rg_off != 0')
L = L.query('rg_off != 0')
tides = tides.query('N != 0')
full = full.query('rg_off != 0')
L['rg_diff_m'] = 2.33*(L['rg_off'] - L['rg_orb'])
R['rg_diff_m'] = 2.33*(R['rg_off'] - R['rg_orb'])
full['rg_diff_m'] = 2.33*(full['rg_off'] - full['rg_orb'])

Rmed = R['rg_diff_m'].rolling('120D').median()
Lmed = L['rg_diff_m'].rolling('120D').median()
(Rmed - Lmed).plot()

# do tides:
import numpy as np
def ENU2rg(E, N, U, heading = 193, inc_angle = 42.27):
    th = np.radians(inc_angle)
    alpha = np.radians(heading)
    LOS = -np.sin(th)*np.cos(alpha)*E + np.sin(th)*np.sin(alpha)*N + np.cos(th)*U
    return LOS


tides['LOS'] = (-1)*ENU2rg(tides['E'], tides['N'], tides['U'])
tides['LOS'] = tides['LOS'] - tides['LOS'][0]
L['rg_diff_m'] = L['rg_diff_m'] - L['rg_diff_m'][0]
L['tidecorrected'] = L['rg_diff_m'] - tides['LOS'] 

# import GACOS data:
import os
import xarray as xr
import rioxarray

from scipy.constants import speed_of_light

def rad2mm_s1(inrad):
    #speed_of_light = 299792458 #m/s
    radar_freq = 5.405e9  #for S1
    wavelength = speed_of_light/radar_freq #meter
    coef_r2m = -wavelength/4/np.pi*1000 #rad -> mm, positive is -LOS
    outmm = inrad*coef_r2m
    return outmm


lat=-20.548
lon=-175.45

epochsdir='/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/102/102D_11050_000002/epochs'

    #wavelength=0.055465763
    #mcoef=4*PI/$wavelength
    #mcoef=226.56
    #gmt grdmath $infile'=gd:Gtiff+n0' 0 NAN $gacos2 $gacos1 SUB 226.56 MUL 0 NAN $U DIV SUB WRAP = $outfile'=gd:Gtiff'

L['gacos'] = 0
gacoscorrs = []
for i,row in L.iterrows():
    gacoscorr = 0
    epoch = int(row['epoch'])
    sltdfile = os.path.join(epochsdir, str(epoch), str(epoch)+'.sltd.geo.tif')
    sltd = rioxarray.open_rasterio(sltdfile)
    #sltd.sel(y=lat, x=lon, method='nearest')
    gacoscorr = float(sltd.median().values)   # this is in radians!
    gacoscorr = rad2mm_s1(gacoscorr)/1000 # in m
    gacoscorrs.append(gacoscorr)


L['gacos'] = gacoscorrs
L['gacos'] = L['gacos'] - L['gacos'][0]
L['gacos_and_tidecorrected'] = L['tidecorrected'] - tides['gacos'] 
    
#def plot_rg(L, name1='rg_diff_m'):
    
fig, ax = plt.subplots()
ax.scatter(L.index, L.rg_diff_m, color='orange', s=10, label='original')
ax.plot(L.index, L.rg_diff_m, color='orange')
ax.scatter(L.index, L.tidecorrected, color='b', s=10, label='SET-corrected')
ax.plot(L.index, L.tidecorrected, color='b')
#ax.plot(L.index, L.tides, color='b', )
ax.set_xlabel('time')
ax.set_ylabel('range offset [m]')
ax.legend()
fig.show()

import urllib
import requests
from lxml import html

def get_ITRF_ENU(lat, lon):
    url = "https://www.unavco.org/software/geodetic-utilities/plate-motion-calculator/plate-motion/model"
    # data to be sent to api
    data = {'name':'modelform',
        'id':'modelform',
        'lat':str(lat),
        'lon':str(lon),
        'model':'itrf2014',
        'format':'ascii'}
    # sending post request and saving response as response object
    r = requests.post(url = url, data = data)
    #outputs are in mm/year, first E, then N
    cont = html.fromstring(r.content)
    [E,N] = cont.text_content().split()[14:16]
    E = float(E)
    N = float(N)
    return E, N

lat=-20.548
lon=-175.45
itrf_e, itrf_n = get_ITRF_ENU(lat, lon)
itrf_los = (-1)*ENU2rg(itrf_e, itrf_n, 0)/1000
timelen = L.index[-1] - L.index[0]
maxpm = itrf_los*timelen.days/365.25

Lmed = L['tidecorrected'].rolling('120D').median()
times = (Lmed.index - Lmed.index[0]).days/365.25
fit = np.polyfit(times,Lmed, 1)
fit_fn = np.poly1d(fit)
fit_value_first = fit_fn(times[0])
fit_value_last = fit_fn(times[-1])
#slope = fit[0] #in m/year

fig, ax = plt.subplots()
#ax.scatter(L.index, L.rg_diff_m, color='orange', s=10, label='original')
#ax.plot(L.index, L.rg_diff_m, color='orange')
ax.scatter(L.index, L.tidecorrected, color='b', s=10, label='SET-corrected rg offsets')
#ax.plot(L.index, L.tidecorrected, color='b')
ax.plot(Lmed.index, Lmed, color='b', label='120 days rolling median')
ax.plot([L.index[0],L.index[-1]], [0, maxpm], color='k', label='ITRF2014 PMM in range: {0:.1f} mm/year'.format(itrf_los*1000))
ax.plot([L.index[0],L.index[-1]], [fit_value_first, fit_value_last], color='g', label='range trend: {0:.1f} mm/year'.format(fit[0]*1000))
#ax.plot(L.index, L.tides, color='b', )
ax.set_xlabel('time')
ax.set_ylabel('range offset [m]')
ax.legend()
fig.show()



cd $slcdir

if [ 1 == 1 ]; then
# and now coregister the slc mosaics (if possible?)
# Co-registration with traditional cross-correlation algorithm
# The procedure for co-registering two SLCs using the cross-correlation algorithm has been
# thoroughly described in the ISP User’s Guide.
# try only classical cross-corr, as cannot use DEM:

for s in `ls`; do if [ ! $m == $s ]; then cd $s;
 spar=$slcdir/$s/$s.$ext.slc.par
 sslc=$slcdir/$s/$s.$ext.slc
 mkdir -p $rslcdir/$s
 rslc=$rslcdir/$s/$s.$ext.rslc
 create_offset $mpar $spar $m'_'$s.off 1 1 1 0
 # init_offset $mslc $sslc $mpar $spar $m'_'$s.off
 init_offset_orbit $mpar $spar $m'_'$s.off
 cp $m'_'$s.off $m'_'$s.off.orbonly
 # offset_pwr $mslc $sslc $mpar $spar $m'_'$s.off $m'_'$s.offs ccp - - offsets 1 - - -
 #offset_pwr $mslc $sslc $mpar $spar $m'_'$s.off $m'_'$s.offs ccp 32 32 offsets 1 64 64 0.1
 offset_pwr $mslc $sslc $mpar $spar $m'_'$s.off $m'_'$s.offs ccp 64 64 offsets 1 16 16 0.1
 #now fit the offsets and improve in off file, if possible.. use only polynomial 1 (translation)!
 offset_fit $m'_'$s.offs ccp $m'_'$s.off coffs coffsets 0.15 1
 # refinement
 offset_pwr $mslc $sslc $mpar $spar $m'_'$s.off $m'_'$s.offs ccp 32 32 offsets 1 32 32 0.1
 offset_fit $m'_'$s.offs ccp $m'_'$s.off coffs coffsets 0.15 1
 # 3b. coregister the mosaics
 # interpolate to RSLC
 SLC_interp $sslc $mpar $spar $m'_'$s.off $rslc $rslc.par - -
cd ..; fi; done
fi



exit

# range pixel offset tracking:
# For a general description of the method
# please refer to Example C in the User’s Guide of the ISP module.
cd $rslcdir
mliwid=`grep range_samples $mpar | gawk {'print $2'}`
mlilen=`grep azimuth_lines $mpar | gawk {'print $2'}`
demwid=226
geolt=/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/102/102D_11050_000002/GEC/gec.lt
dempar=/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/proc/current/102/102D_11050_000002/GEC/fake.seg.dem_par
for s in `ls`; do if [ ! $m == $s ]; then cd $s;
 spar=$rslcdir/$s/$s.$ext.rslc.par
 sslc=$rslcdir/$s/$s.$ext.rslc
 create_offset $mpar $spar tracking.off 1 1 1 0
 # using oversample 2 and then 8x4 (range is 6x higher resolution than azimuth)
 # hmm... no. let's do it 32x8 (i.e. 37x56 m)
 #offset_pwr_tracking $mslc $sslc $mpar $spar tracking.off tracking.offsets tracking.corr 12 4 - 2 -
 offset_pwr_tracking $mslc $sslc $mpar $spar tracking.off tracking.offsets tracking.corr 32 8 - 2 -
 offset_tracking tracking.offsets tracking.corr $mpar tracking.off disp_map disp_val 1 - 0
 widthoff=`grep range_samples tracking.off | awk '{print $2}'`
 lenoff=`grep azimuth_samples tracking.off | awk '{print $2}'`
 # will be only 325 px but let's go with it
 cpx_to_real disp_map disp_map.rng $widthoff 0
 # resample towards orig size
 python3 -c "import cv2; import numpy as np; a = np.fromfile('disp_map.rng', dtype=np.float32).byteswap().reshape(("$lenoff","$widthoff"));cv2.resize(a,dsize=("$mliwid","$mlilen"), interpolation=cv2.INTER_LINEAR).byteswap().tofile('"$s.rng"')"
 #now geocode it ?
 geocode_back $s.rng $mliwid $geolt $s.rng.geo $demwid
 geocode_back disp_map.rng $widthoff $geolt disp_map.rng.geo $demwid
 # and also mli
 multi_look $sslc $spar $s.crop.mli $s.crop.mli.par 1 1
 geocode_back $s.crop.mli $mliwid $geolt $s.crop.mli.geo $demwid
 data2geotiff $dempar $s.rng.geo 2 $s.rng.geo.tif
 data2geotiff $dempar disp_map.rng.geo 2 $s.rng.geo2.tif
 data2geotiff $dempar $s.crop.mli.geo 2 $s.crop.mli.geo.tif
cd ..; fi; done


pip install git+https://github.com/gjoseph92/geogif.git

from geogif import gif
import xarray as xr
... here load all to datacube - don't forget 10*log10
gif(data_array)
quarterly_ndvi = ndvi.resample(time="Q").median()
gif(da, fps=10, to="capecod.gif")




#cosmo:
# do DEM geocoding
  width=`awk '$1 == "range_samples:" {print $2}' ${masterslcdir}/${master}.slc.mli.par`; #echo $width
  length=`awk '$1 == "azimuth_lines:" {print $2}' ${masterslcdir}/${master}.slc.mli.par`; #echo $length
  reducfac=`echo $width | awk '{if(int($1/1000) > 1) print int($1/1000); else print 1}'`
  echo " DEM in radar coordinates does not exist,... computing it "
  demresN=`awk '$1 == "post_lat:" {print $2}' ${dempar}`
  ovrfactN=`echo $demresN $outres | awk '{print -1*($1/$2)}' `
  demresE=`awk '$1 == "post_lon:" {print $2}' ${dempar}`
  ovrfactE=`echo $demresE $outres | awk '{print $1/$2}'`
  echo "  Computing lookup table (from DEM to master geometry) "
gc_map2 $mli.par $dempar $dem ${masterdem}/EQA.dem_par ${masterdem}/EQA.dem ${masterdem}/${master}.lt $ovrfactN $ovrfactE ${masterdem}/ls_map ${masterdem}/ls_map_rdc ${masterdem}/inc ${masterdem}/res ${masterdem}/offnadir ${masterdem}/${master}.sim_sar ${masterdem}/u ${masterdem}/v ${masterdem}/psi ${masterdem}/pix 

width_dem=`awk '$1 == "width:" {print $2}' ${masterdem}/EQA.dem_par`
length_dem=`awk '$1 == "nlines:" {print $2}' ${masterdem}/EQA.dem_par`
  
  echo "  Simulate master amplitude image "
pixel_area ${masterslcdir}/${master}.slc.mli.par ${masterdem}/EQA.dem_par ${masterdem}/EQA.dem ${masterdem}/${master}.lt ${masterdem}/ls_map ${masterdem}/inc ${masterdem}/pix_sigma0 ${masterdem}/pix_gamma0
  #>> $logfile
    echo "  Compute offsets between master amplitude image and simulated master amplitude image "
  create_diff_par ${masterslcdir}/${master}.slc.mli.par - ${masterdem}/${master}.diff_par 1 0 
  #>> $logfile 
init_offset_orbitm $mli ${masterdem}/pix_gamma0 ${masterdem}/$m.diff_par
init_offsetm $mli ${masterdem}/pix_sigma0 ${masterdem}/$m.diff_par
offset_pwr $mslc $sslc $mpar $spar $m'_'$s.off $m'_'$s.offs ccp 64 64 offsets 1 16 16 0.1

offset_pwrm ${masterdem}/pix_sigma0 ${masterslcdir}/${master}.slc.mli ${masterdem}/${master}.diff_par ${masterdem}/${master}.offs ${masterdem}/${master}.ccp 256 256 offsets 2 64 64 0.2
  # >> $logfile
  offset_fitm ${masterdem}/${master}.offs ${masterdem}/${master}.ccp ${masterdem}/${master}.diff_par ${masterdem}/coffs ${masterdem}/coffsets 0.2 1 
  # >> $logfile
  echo "  Refine lookup table using offsets "
  gc_map_fine ${masterdem}/$master.lt ${width_dem} ${masterdem}/$master.diff_par ${masterdem}/$master.lt_fine 1  >> $logfile
  echo "  Geocode the master amplitude image "
  geocode_back ${masterslcdir}/${master}.slc.mli $width ${masterdem}/$master.lt_fine ${masterdem}/EQA.${master}.slc.mli ${width_dem} ${length_dem} 2 0  >> $logfile
  echo "  Convert to master radar geometry the cropped DEM (geo/EQA.dem): geo/${master}.hgt "
  
  
  geocode ${masterdem}/$master.lt_fine ${masterdem}/EQA.dem ${width_dem} ${masterdem}/${master}.hgt ${width} ${length} 2 0 >> $logfile
  
  


geocode_back $mli $mliwid lt mli.geo $demwid    # to geocode MLI
geocode ${masterdem}/$master.lt_fine ${masterdem}/EQA.dem ${width_dem} ${masterdem}/${master}.hgt ${width} ${length} 2 0 

# do coreg - done already, no problemo

#make ifgs:
for x in `ls SLC`; do
mkdir IFG/$x; create_diff_par SLC/20190603/20190603.slc.par RSLC/$x/$x.rslc.par IFG/$x/diffpar 0 0;
phase_sim_orb SLC/20190603/20190603.slc.par RSLC/$x/$x.rslc.par IFG/$x/diffpar - IFG/$x/simpha SLC/20190603/20190603.slc.par
SLC_diff_intf SLC/20190603/20190603.slc RSLC/$x/$x.rslc SLC/20190603/20190603.slc.par RSLC/$x/$x.rslc.par IFG/$x/diffpar IFG/$x/simpha IFG/$x/diff 1 1
cpxfiddle -w 16286 -f cr4 -q phase -M 16/16 -o sunraster -c jet -Bb IFG/$x/diff > IFG/$x/diff.ras
convert IFG/$x/diff.ras IFG/$x/diff.png; rm IFG/$x/diff.ras


create_diff_par $mlipar geo/dem_seg_par geo_crop/$m.diff_par 1 0
init_offsetm $mli geo/sim_sar geo_crop/$m.diff_par
init_offset_orbitm $mli geo/sim_sar geo_crop/$m.diff_par
create_diff_par SLC/20190603/20190603.slc.par RSLC/$x/$x.rslc.par IFG/$x/diffpar 1 0
phase_sim_orb SLC/20190603/20190603.slc.par RSLC/$x/$x.rslc.par IFG/$x/diffpar - IFG/$x/simpha SLC/20190603/20190603.slc.par

done



'
fucking hell with gamma. do cosmooo! fuck
