#!/bin/bash
###########################################################
# For CEMS configuration
###########################################################
# #!/bin/bash
 grouppath=/group_workspaces/cems2/nceo_geohazards
 LiCSARpath=${grouppath}/software/LiCSAR2.0
 # External information
 export DEMs_DIR=${grouppath}/DEMs   # For the moment only works for srtm_1arcsec with filled voids
 export ORBs_DIR=${grouppath}/orbits/s1a
 # Load the bash library
 source ${LiCSARpath}/lib/LiCSAR_bash_lib.sh     # Library of bash functions
 # Load GAMMA software
 export PATH="${grouppath}/software/GAMMA/20151209/2/default/bin.exec:${PATH}"
 # Update the PATH to contain main programs
 export PATH=${LiCSARpath}/bin/:$PATH      # Main scripts to command flow
 export PATH=${LiCSARpath}/scripts/:$PATH  # Secondary scripts
 export PATH=${LiCSARpath}/external/:$PATH # Third party programs (e.g, fastjar if jar is not available)
 export PATH=${LiCSARpath}/python/:$PATH   # python programs: burst-related, pyaps interface scripts,... etc.
 # Add pyAPS package to PYTHONPATH to use as a library [atmospheric corrections]
 export PYTHONPATH=${LiCSARpath}/python:${LiCSARpath}/lib:$PYTHONPATH
###########################################################

# [basic User Params]
#export archpath = /nfs/a1/raw/sentinel/
#export S1Aorbitpath = /nfs/a1/raw/orbits/S1A/

# [sql server]
host = foe-db.leeds.ac.uk
dbname = sees1ice
#dbname = sees1archive
dbuser = s1user
dbpass = Aerahf7cee2j

# [debug]
# If there are issues, use this to provide extra output
# Level 0 is least verbose
# Level 1 will output information about current processes
# Level 2 will output information about variables and settings
# Level 3 will include the verbose HTTP connection info, in files found in ../tmpfiles.  
# These files need manual deletion once debugging complete
# Level 4 can be used to print output rather than download (in which case dl = 3 won't return anything)
dl = 2

# [lockfile]
# Define the details for a hidden lockfile to be created at the start of a download and deleted at the end
# Create the path on a networked server to ensure that only one instant is running at once
# DO NOT ALTER this unless you know what you are doing
lockPath = /nfs/a147/insar/sentinel1/S1_archive/
lockFile = .S1_DandA.lock
