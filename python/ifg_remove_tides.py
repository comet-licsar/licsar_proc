import sys
import rioxarray as r
import numpy as np

hgtfile=sys.argv[1]
ifile=sys.argv[2]
t1f=sys.argv[3]
t2f=sys.argv[4]
ofile=sys.argv[5]

h=r.open_rasterio(hgtfile)
i=r.open_rasterio(ifile)
t1=r.open_rasterio(t1f)
t2=r.open_rasterio(t2f)

h.values=i.values
dt=(t1-t2)*226.56
h=h-dt
h.values=np.angle(np.exp(1j*h.values))

h.rio.to_raster(ofile)
