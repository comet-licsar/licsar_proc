# this is set of functions developed to get nice AHB paper figures
import pygmt
import os
import pandas as pd
import numpy as np
import datetime as dt
import xarray as xr


fr='005D_05199_131313'

ahbdir='/gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/AHB'
mahb=os.path.join(ahbdir, 'Milan')
figdir=os.path.join(ahbdir+'_figures', 'ionoetc')
outfile=os.path.join(ahbdir,'ionograds',fr+'.nc')


# now this will perform iono gradient estimation



h5file = os.path.join(ahbdir, fr, fr+'.cum.h5')
a=xr.open_dataset(h5file) #, engine="h5netcdf")


g='iono'
def gradts(a,g,resolm=300):
    vals = []
    for f in a[g]:
        gradd=np.gradient(f.values)
        gradd=np.sqrt(gradd[0]**2 + gradd[1]**2)
        vals.append(np.nanmean(gradd))
    vals = np.array(vals)
    return vals / resolm

resoldeg = float(a.post_lon.values)
resolm=resoldeg*111111

# first refer to 2020-02-01 as approx in the centre of (all) the dataset(s) + very low sunspot no.
imd=list(a.imdates.values)+[20200201]
imd=np.array(imd)
i=np.argmax(np.argsort(imd))-1
a['iono']=a['iono']-a['iono'][i]

# and actually drop that reference as it is zero (would degrade the mean?)
# or actually not, let's leave it...
# iono = a['iono']
#iono = np.concatenate([iono[:i], iono[i+1:]], axis=0)
# imdates = a['imdates']
#imdates = np.concatenate([imdates[:i], imdates[i+1:]], axis=0)

# now get the iono gradients
valsI = gradts(a,'iono', resolm)
#valsT = gradts(a,'tide')

x=a.imdates.values
xx=[]
for xa in x:
    xx.append(pd.Timestamp(str(xa)))


x = np.array(xx)

# do monthly avg
df = pd.DataFrame(
    data={
        'x': x[1:],
        'ionograd': valsI[1:]*1000
    }
)
# Convert the date column to a datetime object if it's not already
df["date"] = pd.to_datetime(df["x"])

# Group by year and month, and calculate mean and std
monthly_stats = df.groupby(df["date"].dt.to_period("M"))["ionograd"].agg(["mean", "std"])
sss=monthly_stats['std'].values
sss[np.isnan(sss)]=0
monthly_stats['std'] = sss

# Convert PeriodIndex back to datetime (optional, for better visualization)
monthly_stats.index = monthly_stats.index.to_timestamp()


# now add it to the final xr.Dataset?
# or just... export to csv..
meanlat = float(a.corner_lat+a.post_lat*a.phony_dim_2.mean())
months = pd.date_range(start="2016-03", periods=12 * 8 + 1, freq="MS")  # Monthly start dates

monthly_stats = monthly_stats.reindex(months, fill_value=np.nan)

dataset = xr.Dataset(
    {
        "meangrad": (["month"], monthly_stats['mean']),
        "stdgrad": (["month"], monthly_stats['std'])
    },
    coords={
        "month": months,
        #"meanlat": meanlat
        #"meanlat_stdgrad": meanlat_stdgrad
    }
)
dataset.attrs['unit']='mm/km'
dataset.attrs['meanlat']=meanlat
dataset.attrs['orbit']=fr[3]
dataset.attrs['frame']=fr


dataset.to_netcdf(outfile)