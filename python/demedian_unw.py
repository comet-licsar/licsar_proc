#!/usr/bin/env python
import numpy as np
import sys, os, fnmatch

filename=sys.argv[1]
from LiCSAR_lib.unwrp_lib import demedian_unw

demedian_unw(filename)
