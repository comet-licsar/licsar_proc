#!/usr/bin/env python

# this will output volcano name from the volclip

import volcdb as volc
import sys

volclipid=int(sys.argv[1])

vname = volc.get_volc_info(volc.get_volcano_from_vid(volclipid)).name.values[0].replace(' ','_')
print(vname)
