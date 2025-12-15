#!/usr/bin/env python

# first of all, the nc files must be generated using ahb_extract_tide_grad.py and ahb_extract_iono_grad.py


'''
create TS of iono corrections:
sunspots (monthly average with std?)
TEC per used frame/epoch/time, coloured by latitude
correction max gradient
Materials for this figure:
[Data]
https://www.sidc.be/SILSO/datafiles#total

[pygmt]
https://www.pygmt.org/latest/gallery/lines/envelope.html#sphx-glr-gallery-lines-envelope-py
https://www.pygmt.org/latest/gallery/lines/line_custom_cpt.html#sphx-glr-gallery-lines-line-custom-cpt-py
'''


import pygmt
import os
import pandas as pd
import numpy as np
import datetime as dt
import xarray as xr
#from lics_unwrap import *


ahbdir='/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/AHB'
mahb=os.path.join(ahbdir, 'Milan')
figdir=os.path.join(ahbdir+'_figures', 'ionoetc')
#os.listdir(figdir)


# 1. prepare sunspot data
monthlycsv = os.path.join(figdir, 'SN_m_tot_V2.0.csv')
sunspots = pd.read_csv(monthlycsv, sep=';', header=None)
#Column 3: Date in fraction of year.
#Column 4: Monthly mean total sunspot number.
#Column 5: Monthly mean standard deviation of the input sunspot numbers.

sunspots = sunspots[sunspots[0]>=2014]
sunspots = sunspots[sunspots[0]<2025]

xsun = sunspots[2].values
ysun = sunspots[3].values
stdsun = sunspots[4].values

xsundt = []
from datetime import datetime, timedelta
for year_decimal in xsun:
    # Extract the year and the fractional part
    year = int(year_decimal)
    fraction = year_decimal - year
    # Calculate the number of days corresponding to the fraction
    days_in_year = 365 + (1 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 0)  # Check for leap year
    days = fraction * days_in_year
    # Create the datetime object
    date = (datetime(year, 1, 1) + timedelta(days=days)).date()
    xsundt.append(date)

# Define a pandas.DataFrame with columns for x and y as well as the lower and upper
# deviations
df_sun = pd.DataFrame(
    data={
        'x': xsundt,
        'y': ysun,
        'y_deviation_low':stdsun,
        'y_deviation_upp':stdsun
    }
)


# for iono :


import glob
descs = glob.glob('/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/AHB/ionograds/*D*.nc')
ascs = glob.glob('/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/AHB/ionograds/*A*.nc')

i=0
ncpath = descs[i]
def get_x_y_lat(ncpath):
    gradxr=xr.open_dataset(ncpath)
    x = gradxr.month.values
    y = gradxr.meangrad.values
    lat = gradxr.attrs['meanlat']
    return x,y,lat


fig = pygmt.Figure()
pygmt.config(FORMAT_DATE_MAP="yyyy")
pygmt.config(FORMAT_TIME_MAP="yyyy")
pygmt.config(FORMAT_TIME_PRIMARY_MAP="yyyy")

projection="X18c/6c"
#region=[dt.datetime(2014,10,1).date(), dt.datetime(2024,12,31).date(), 0, np.max(ysun)+np.max(stdsun)]
region=[dt.datetime(2015,1,1).date(), dt.datetime(2024,12,31).date(), 0, np.max(ysun)+np.max(stdsun)]
regiongrad=[dt.datetime(2015,1,1).date(), dt.datetime(2024,12,31).date(), 0, 0.7]
# Specific date to add the vertical dashed line
vertical_date = pd.Timestamp("2020-02-01")

import glob
descs = glob.glob('/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/AHB/ionograds/*D*.nc')
ascs = glob.glob('/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/AHB/ionograds/*A*.nc')

#####################################
fig.basemap(
    region=region,
    #projection="X10c",
    projection=projection,
    frame=["SEn", "xa1yf1y+lyear", "ya100f50+lsunspot number"] #, "+t\"%Y\""]
    #frame=["WSne+tsymmetric deviations +d"] #, "xa2f1", "ya1f0.1"],
)


fig.plot(
    data=df_sun,
    close="+d",
    # Fill the envelope in gray color with a transparency of 50 %
    fill="gray@50",
    pen="1p,gray30",
    label="monthly average sunspot number"
)

fig.basemap(
    region=regiongrad,
    #projection="X10c",
    projection=projection,
    frame=["W", "yaf+liono gradient [mm/km]"] #, "+t\"%Y\""]
    #frame=["WSne+tsymmetric deviations +d"] #, "xa2f1", "ya1f0.1"],
)

gradxr=xr.open_dataset(descs[i])
gradxr=xr.open_dataset(ascs[i])
df_grad_des =  pd.DataFrame(
    data={
        'x': gradxr.month.values,
        'y': gradxr.meangrad.values,
        'y_deviation_low':gradxr.stdgrad.values,
        'y_deviation_upp':gradxr.stdgrad.values
    }
)
fig.plot(
    data=df_grad_des,
    #close="+d",
    # Fill the envelope in gray color with a transparency of 50 %
    #fill="red@70",
    #pen="1p,red",
    pen="thick,green",
    label="mean iono correction gradient" # in descending frames"
)

pygmt.makecpt(cmap="roma", series=[20, 50])
for i in range(10): # range(len(descs)):
    #x, y, lat = get_x_y_lat(descs[i])
    x, y, lat = get_x_y_lat(ascs[i])
    fig.plot(x=x, y=y, cmap=True, zvalue=lat, pen="thick,+z,-")

# Add a vertical dashed line
fig.plot(
    x=[vertical_date, vertical_date],  # X positions (start and end of the line are the same for verticality)
    y=[0, regiongrad[-1]/2],  # Y positions (full height of the plot)
    pen="1p,black,--"  # Line style: thickness 1p, color black, dashed (--)
)

fig.legend(position="JTL+jTL+o0.2c", box=False)
fig.colorbar(position="JTL+jTL+o1.1c+h+w4c", frame=["xaf+llatitude"])



















############ SAME BUT FOR GACOS
fig = pygmt.Figure()
pygmt.config(FORMAT_DATE_MAP="yyyy")
pygmt.config(FORMAT_TIME_MAP="yyyy")
pygmt.config(FORMAT_TIME_PRIMARY_MAP="yyyy")

projection="X9c/5c"
#region=[dt.datetime(2014,10,1).date(), dt.datetime(2024,12,31).date(), 0, np.max(ysun)+np.max(stdsun)]
#region=[dt.datetime(2015,1,1).date(), dt.datetime(2024,12,31).date(), 0, np.max(ysun)+np.max(stdsun)]
regiongrad=[dt.datetime(2015,1,1).date(), dt.datetime(2024,12,31).date(), 0, 12]
regiongrad=[dt.datetime(2016,3,1).date(), dt.datetime(2024,3,1).date(), 0, 12]
region=regiongrad
# Specific date to add the vertical dashed line
vertical_date = pd.Timestamp("2020-02-01")

import glob
descs = glob.glob('/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/AHB/gacosgrads/*D*.nc')
ascs = glob.glob('/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/AHB/gacosgrads/*A*.nc')

doscs = descs
doscs = ascs
'''
#####################################
fig.basemap(
    region=region,
    #projection="X10c",
    projection=projection,
    frame=["SEn", "xa1yf1y+lyear", "ya100f50+lsunspot number"] #, "+t\"%Y\""]
    #frame=["WSne+tsymmetric deviations +d"] #, "xa2f1", "ya1f0.1"],
)


fig.plot(
    data=df_sun,
    close="+d",
    # Fill the envelope in gray color with a transparency of 50 %
    fill="gray@50",
    pen="1p,gray30",
    label="monthly average sunspot number"
)
'''
fig.basemap(
    region=regiongrad,
    #projection="X10c",
    projection=projection,
    frame=["SeWn", "xa1yf1y+lyear", "yaf+ltropo gradient [mm/km]"] #, "+t\"%Y\""]
    #frame=["WSne+tsymmetric deviations +d"] #, "xa2f1", "ya1f0.1"],
)

i=0
gradxr=xr.open_dataset(descs[i])
gradxr=xr.open_dataset(ascs[i])
df_grad_des =  pd.DataFrame(
    data={
        'x': gradxr.month.values,
        'y': gradxr.meangrad.values,
        'y_deviation_low':gradxr.stdgrad.values,
        'y_deviation_upp':gradxr.stdgrad.values
    }
)
fig.plot(
    data=df_grad_des,
    #close="+d",
    # Fill the envelope in gray color with a transparency of 50 %
    #fill="red@70",
    #pen="1p,red",
    pen="thick,green",
    label="mean tropo correction gradient" # in descending frames"
)

pygmt.makecpt(cmap="roma", series=[20, 50])
for i in range(len(doscs)):
    #x, y, lat = get_x_y_lat(descs[i])
    x, y, lat = get_x_y_lat(doscs[i])
    fig.plot(x=x, y=y, cmap=True, zvalue=lat, pen="thin,+z,-")

# Add a vertical dashed line
fig.plot(
    x=[vertical_date, vertical_date],  # X positions (start and end of the line are the same for verticality)
    y=[0, regiongrad[-1]/2],  # Y positions (full height of the plot)
    pen="1p,black,--"  # Line style: thickness 1p, color black, dashed (--)
)

#fig.legend(position="JTL+jTL+o0.2c", box=False)
#fig.colorbar(position="JTL+jTL+o1.1c+h+w4c", frame=["xaf+llatitude"])

# Plot the data points on top
#fig.plot(data=df_sun, style="c0.1c", pen="1p,gray30", fill="darkgray")

fig.show()