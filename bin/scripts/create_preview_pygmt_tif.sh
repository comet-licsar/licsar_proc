#!/bin/bash
fname=$1
python3 -c "import lics_processing as lp; import lics_vis as lv; gg=lv.pygmt_plot(lp.load_tif2xr('"$fname"'), title='"$fname"', lims=None, projection='M10c', label='units'); gg.savefig('"$fname".png', dpi=120)" 2>/dev/null
ls $fname".png"
