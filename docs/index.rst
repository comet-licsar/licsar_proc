LiCSAR proc
===========

This section describes core functionality of LiCSAR system, that is used by LiCSAR FrameBatch system but can be used independently.
Some parameters of particular use (often not fully integrated in FrameBatch) are discussed here, yet we assume that you work with
standard project started by licsar_make_frame.sh.

Forming LiCSAR interferograms (Sentinel-1)
------------------------------

The core processing scripts are also described (in more details, and still valid, although of older date) in https://gitlab.com/comet_licsar/licsar_documentation/-/wikis/ciw2019/licsar

1. create SLC files
^^^^^^^^^^^^^^^^^^^^^^^

see
::
  LiCSAR_01_mk_images.py -h

normal use:
::
  frameid= # e.g. 016A_02562_091313
  m= # reference date, e.g. 20191103
  s= # e.g. 20220101
  e= # e.g. 20220201
  LiCSAR_01_mk_images.py -n -d . -f $frameid -m $m -s $s -e $e -a 4 -r 20

This command would create SLC for given frame in given time period. For this to work, you should be sure the data are existing on disk and are
ingested to LiCSInfo database (command arch2DB.py). This is done automatically through licsar_make_frame.sh.

2. coregister to RSLC files
^^^^^^^^^^^^^^^^^^^^^^^

see:
::
  LiCSAR_02_coreg.py -h

normal use:
::
  frameid= # e.g. 016A_02562_091313
  m= # reference date, e.g. 20191103
  LiCSAR_02_coreg.py -f $frameid -d . -m $m -i

This command would coregister/resample all files in SLC folder to RSLC directory, for given frame.
The standard constrains apply, for example limit for spectral diversity estimation to combine epochs with up to Btemp=180 days.


In order to force-skip the 180 days limit, use parameter '-E' (some extra tweaks there, EIDP-related, including skip of Btemp check).
In practice, to connect a cluster of epochs far in time, choose one SLC in that cluster that has no missing bursts (see e.g. their mli preview, or check size),
and add its date to slclist.txt and try coregister that one by ``LiCSAR_02_coreg.py -f $frameid -d . -m $m -i -l slclist.txt -E``.


Only then start the standard procedure (without -E) for all the SLCs left. This way you would make sure that the other epochs would use appropriate
RSLC3 for the spectral diversity (see e.g. https://www.mdpi.com/2072-4292/12/15/2430 )


3. create interferograms
^^^^^^^^^^^^^^^^^^^^^^^

see:
::
  LiCSAR_03_mk_ifgs.py -h

The script is very useful if you have your own list of interferograms to form, e.g. in a text file containing lines as '20200101_20200202' etc. (see help).
An extra parameter -n would use parallel processing on given number of CPUs.


4. unwrap interferograms
^^^^^^^^^^^^^^^^^^^^^^^

For the original unwrapping approach, running on radar-coordinate interferograms, use:
::
  LiCSAR_04_unwrap.py -h

and then you may geocode the result, as discussed in next section.

However, you may find useful (and faster) the updated version, currently used by LiCSAR FrameBatch, that performs unwrapping on already geocoded wrapped interferograms:
::
  unwrap_geo.sh

Finally, you may experiment with the updated (much improved) unwrapper, running through python, and starting again from geocoded interferograms. This script is used by licsar2licbas.sh described later.
::
  import unwrp_multiscale as unw
  help(unw.process_frame)


5. geocoding results
^^^^^^^^^^^^^^^^^^^^^^^
For geocoding results, please use the following command:
::
  create_geoctiffs_to_pub.sh


Post-processing
-------------------

LiCSAR to LiCSBAS (JASMIN)
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
  #-s ....... if the reunwrapping is to be performed, use smoothing (two-pass unw approach, similar effect as with cascade, only milder)
  #-H ....... this will use hgt to support unwrapping (only if using reunwrapping)
  #-T ....... use testing version of LiCSBAS
  #-S ....... strict mode - e.g. in case of GACOS, use it only if available for ALL ifgs
  #-G lon1/lon2/lat1/lat2  .... clip to this AOI
  ##
  ## following is an ongoing work, for testing only:
  ##-C ....... use coherence stability index instead of orig coh per ifg (experimental - might help against loop closure errors, maybe)
  ##-k ....... use cohratio everywhere (i.e. for unwrapping, rather than orig coh - this is experimental attempt)



While parameters -C, -k are only related to a short-term experiment (should conclude in use of amplitude stability and/or general coherence for masking and weighting),
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


Decomposition to E-U(+N) vectors
^^^^^^^^^^^^^^^^^^^^^^^

This section should contain information on both decomposition from A+D (use of Andrew's tutorial?)

Bringing ENU model values to line-of-sight
^^^^^^^^^^^^^^^^^^^^^^^

Inverse procedure (with example?) using E,N,U tif files to convert ENU->LOS.
