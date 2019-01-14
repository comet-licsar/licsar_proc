#!/bin/bash
###########################################################
# [basic User Params]
#export archpath = /nfs/a1/raw/sentinel/
#export S1Aorbitpath = /nfs/a1/raw/orbits/S1A/

# [sql server]
# this is for CIPEG db that is actually not used anymore:
##host = foe-db.leeds.ac.uk
#host = see-cipegdb.leeds.ac.uk
#dbname = sees1ice
##dbname = sees1archive
#dbuser = s1user
#dbpass = Aerahf7cee2j
# functional db in CEMS for COMET LiCS db
host = 130.246.129.18
dbname = licsinfo_live
dbuser = lics
dbpass = T34mLiCS

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
