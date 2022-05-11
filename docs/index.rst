LiCSAR proc
===========

This section describes core functionality of LiCSAR system, that is used by LiCSAR FrameBatch system but can be used independently.
Some parameters of particular use (often not fully integrated in FrameBatch) are discussed here.

1. create SLC files
-------------------

test of the integrated docs

2. coregister to RSLC files
-------------------

how to coreg outside of 180 days limit? .... just run LiCSAR_02_coreg.py with '-E' parameter (some extra tweaks there, EIDP-related)

test of the integrated docs

3. create interferograms
-------------------

test of the integrated docs

4. unwrap interferograms
-------------------

test of the integrated docs

5. additional functionality
-------------------

LiCSAR to LiCSBAS
^^^^^^^^^^^^^^^^^^^^^^^
This script runs LiCSBAS processing from the LiCSAR data. To be used in JASMIN environment.

The script would read frame data from $LiCSAR_public directory, prepare them for LiCSBAS and run LiCSBAS with default parameters.
If you run the script from directory with your GEOC outputs, it would instead use the local data from this folder.
Afterwards, you may just fine tune parameters of LiCSBAS step 15 (and 16) and rerun them, for the final result.
::
  licsar2licsbas.sh frame [startdate] [enddate]
  #e.g. 155D_02611_050400 20141001 20200205
  #parameters:
  #-M 10 .... this will do extra multilooking (in this example, 10x multilooking)
  #-u ....... use the (extra Gaussian-improved multilooking and) reunwrapping procedure (useful if multilooking..)
  #-c ....... if the reunwrapping is to be performed, use cascade (might be better, especially when with shores)
  #-s ....... use coherence stability index instead of orig coh per ifg (experimental - might help against loop closure errors, maybe)
  #-k ....... use cohratio everywhere (i.e. for unwrapping, rather than orig coh - this is experimental attempt)
  #-H ....... this will use hgt to support unwrapping (only if using reunwrapping)
  #-T ....... use testing version of LiCSBAS
  #-S ....... strict mode - e.g. in case of GACOS, use it only if available for ALL ifgs
  #-G lon1/lon2/lat1/lat2  .... clip to this AOI


While parameters -s, -k are only related to a short-term experiment (should conclude in use of amplitude stability and/or general coherence for masking and weighting),
the other parameters are practically used/recommended to understand.


Explaining on example, use of
::
  licsar2licsbas.sh -c -M 5 -u -T -G 5.1/5.2/3.3/3.5 100D_00000_010101 20150101 20160101

would grab **wrapped** interferograms of this (fictive) frame 100D that cover period of year 2015, then it will check for availability of GACOS corrections and use them if they exist for most of epochs
(if you used -S, GACOS corrections would be applied only if they exist for ALL epochs). Then it would crop them to the coordinates given by -G, and then it will **reunwrap** them (-u) with 5x multilooking
(so the resolution if using default LiCSAR data would become approx. 500 m), with support of cascade approach (-c) that means a longer wave signal is first estimated/unwrapped (using 10x the -M factor)
and used to bind the final unwrapped result - therefore especially decorrelated areas would not induce unwrapping error.. hopefully.

The data here will be prepared to folder GEOCml5GACOSclip.
Then, the -T would use up-to-date LiCSBAS codes with their experimental functionality ON (in this case, e.g. nullification of pixels in unwrapped pairs with loop closure errors over pi is ON).
The whole procedure will run in the background through JASMIN's LOTUS server (see generated .sh files) and once finished, results will be in TS_GEOCml5GACOSclip, plus additional files will be generated
(e.g. geotiffs of velocity estimate, or standard NetCDF file that can be loaded to e.g. QGIS or ncview to plot time series from 'cum' layer, etc.)
