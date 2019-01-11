import os
import sys
import subprocess as subp

config = {
    'VERSION': 'V2.0',
    'ENV': 'development', #[production|development|mirror]
    'DEST': 'leeds' #[CEMS|leeds]
}

rglks = 20
azlks = 4
outres = 0.001
batchflag=False

configfile = os.environ["LiCSARconfig"]

# if config['DEST'] == 'CEMS':
    # # Due to file size check on CEMS and the hardware reporting inaccurate available disk space, a temp file write is needed
    # os.environ["LiCSAR_temp"] = "/work/scratch/" + os.environ["USER"]
    # if not os.path.exists(os.environ["LiCSAR_temp"]):
        # os.makedirs(os.environ["LiCSAR_temp"])

    # # # appends the gamma software bin path to the current process encironemt path, will clear at program end.
    # # os.environ["PATH"] += os.pathsep + '/group_workspaces/cems2/nceo_geohazards/software/GAMMA/20151209/bin'
    # os.environ["grouppath"] = "/group_workspaces/cems2/nceo_geohazards"
    # os.environ["LiCSARpath"] = LiCSAR_PATH

    # # External information
    # os.environ["DEMs_DIR"] = os.environ["grouppath"] + "/DEMs"   # For the moment only works for srtm_1arcsec with filled voids
    # os.environ["ORBs_DIR"] = os.environ["grouppath"] + "/orbits/s1a"

    # # Load the bash library
    # # Need to add to the sys environment the source of the licsar bash scripts
    # output = subp.check_output("source " + LiCSAR_PATH + "/lib/LiCSAR_bash_lib.sh; env -0", shell=True,
                               # executable="/bin/bash")
    # os.environ.update(l.partition('=')[::2] for l in output.split('\0'))

    # # Load GAMMA software
    # os.environ["PATH"] = os.environ["grouppath"] + "/software/GAMMA/20151209/2/default/bin.exec" + os.pathsep + os.environ["PATH"]

    # # Load snaphu software
    # os.environ["PATH"] = os.environ["grouppath"] + "/software/snaphu/bin" + os.pathsep + os.environ["PATH"]

    # # Add pyAPS package to PYTHONPATH to use as a library [atmospheric corrections]
    # if "PYTHONPATH" in os.environ.keys():
        # os.environ["PYTHONPATH"] = LiCSAR_PATH + "/python" + os.pathsep + LiCSAR_PATH +"/lib" + os.pathsep + os.environ["PYTHONPATH"]
    # else:
        # os.environ["PYTHONPATH"] = LiCSAR_PATH + "/python" + os.pathsep + LiCSAR_PATH + "/lib"


