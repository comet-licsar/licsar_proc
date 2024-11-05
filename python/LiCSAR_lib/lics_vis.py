import pygmt
import matplotlib.pyplot as plt
import numpy as np
try:
    import contextily as ctx
except:
    print('contextily is not installed')

def plot3(A,B,C):
    '''inputs are three xr.dataarrays to plot'''
    origfigsize = plt.rcParams['figure.figsize']
    plt.rcParams["figure.figsize"] = [18,4]
    plt.subplot(1,3,1)
    #AA.origpha
    A.rename('rad').plot()
    plt.subplot(1,3,2)
    #AA.unwlow
    B.rename('rad').plot()
    plt.subplot(1,3,3)
    #AA.toremove
    C.rename('rad').plot()
    plt.show()
    plt.rcParams['figure.figsize']=origfigsize


def licsbas_pygmt_plot(cube, title = 'vel', vminmax = [-15, 15],
                       toplot = 'vel', volcid = None):
    ''' plots licsbas result in netcdf, loaded to xr.Dataset (cube).
    if volcano ID is set, it will add its symbol to the plot.
    '''
    import geopandas as gpd
    from shapely.geometry import Point
    refpoint = Point((cube.ref_lon, cube.ref_lat))
    refpoint = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[refpoint])

    #print('Plotting velocity layer on DEM')
    grid = cube[toplot]
    if maskit:
        grid = grid.where(cube.mask == 1)
    if volcid:
        import volcdb as v
        vid = v.get_volclip_vids(volcid)[0]
        vgpd = v.get_volclips_gpd(vid)
        # vgpd.geom.bounds
        minlon = float(vgpd.geom.bounds.minx)
        minlat = float(vgpd.geom.bounds.miny)
        maxlon = float(vgpd.geom.bounds.maxx)
        maxlat = float(vgpd.geom.bounds.maxy)
        grid = grid.sel(lon=slice(minlon, maxlon), lat=slice(minlat, maxlat))
    #grid = grid.load()
    gridgmt = pygmt_plot(grid, title, 'mm/year', vminmax)
    gridgmt.plot(data=refpoint, pen="4p,magenta4")
    if volcid:
        volctb = v.get_volc_info(volcid)
        gridgmt.plot(data=volctb.geom, style="t8p", pen="0.5p,black", fill="red")



def volcano_clip_plot(volcid, bevel = 0.1):
    '''plots volcano given by its volcano id - see volcdb'''
    #volcid = 243040
    import volcdb as volc
    volcrecord = volc.get_volc_info(volcid)
    vid = volc.get_volclip_vids(volcid)[0]
    volclip = volc.get_volclips_gpd(vid)
    # CTX
    sourcetiles = ctx.providers.Esri.WorldImagery
    lon1, lon2, lat1, lat2 = float(volclip.bounds.minx), float(volclip.bounds.maxx), float(volclip.bounds.miny), float(
        volclip.bounds.maxy)
    fig = pygmt.Figure()
    region = [lon1 - bevel, lon2 + bevel, lat1 - bevel, lat2 + bevel]
    fig.tilemap(
        region=region,  # projection=projection,
        # region=[-157.84, -157.8, 21.255, 21.285],
        # projection="M12c",
        # Set level of details (0-22)
        # Higher levels mean a zoom level closer to the Earth's
        # surface with more tiles covering a smaller
        # geographic area and thus more details and vice versa
        # Please note, not all zoom levels are always available
        zoom=14,
        # Use tiles from OpenStreetMap tile server
        source=sourcetiles
    )
    fig.coast(region=region, shorelines=True)  # , water="lightblue")
    fig.plot(volclip.geom)
    fig.plot(volcrecord.geom, color='red')
    # fig.plot(gpd_overlaps)
    return fig #fig.show()


def pygmt_plot_interactive(grid, title, label='deformation rate [mm/year]', lims=[-25, 10],
                           cmap="roma", photobg=False, plotvec=None):
    import IntPyGMT.IntPyGMT_overlay as ingmt
    # print("%matplotlib widget")
    tempng = '/tmp/pygmt_pi.png'
    fig, region, projection, xshift, yshift = pygmt_plot(grid, title, label, lims,
               cmap, photobg, plotvec, interactive=True)
    fig.savefig(tempng)
    ingmt.gmt_png(tempng, region, projection, xshift, yshift)


def pygmt_plot(grid, title, label='deformation rate [mm/year]', lims=[-25, 10],
               cmap="roma", photobg=False, plotvec=None, interactive = False):
    ''' Function to generate (nice) plot of given grid using pyGMT
    
    Args:
        grid (xr.DataArray): input grid to plot
        title (str):  title (note too long title will disbalance the figure)
        label (str):  label below the colour scale
        lims (list):  colour scale limits
        cmap (str):   colour scale map (try 'vik' for E-W)
        photobg (bool): will plot orthophotomap as the background (if False, DEM relief is used)
        plotvec (geopandas etc): will plot vector data to the map, using pyGMT defaults
    
    Returns:
        pygmt.figure.Figure
    '''
    try:
        grid = grid.load()
    except:
        print('error loading the input dataarray to memory')
        return False
    # try cmap 'vik' for E-W
    #
    # grid = a['U'].where(a.mask < 5) - 10
    # topo_data = '@earth_relief_03s' #3 arc second global relief (SRTM3S)
    topo_data = '@earth_relief_01s'  # 3 arc second global relief (SRTM3S)

    minlon, maxlon = float(np.min(grid.lon)), float(np.max(grid.lon))
    minlat, maxlat = float(np.min(grid.lat)), float(np.max(grid.lat))


    fig = pygmt.Figure()
    pygmt.config(FORMAT_GEO_MAP="ddd.xx") #, MAP_FRAME_TYPE="plain")
    projection = "M13c" # 'R13c' for Robinson etc.
    region = [minlon, maxlon, minlat, maxlat]

    if interactive:
        xshift = '1.5c'
        yshift = '2.5c'
        fig.shift_origin(xshift=xshift, yshift=yshift)
    fig.basemap(region=region, projection=projection, frame=["af", '+t"{0}"'.format(title)])

    if photobg:
        import contextily as ctx
        sourcetiles = ctx.providers.Esri.WorldImagery
        fig.tilemap(
            region=region, projection=projection,
            # region=[-157.84, -157.8, 21.255, 21.285],
            # projection="M12c",
            # Set level of details (0-22)
            # Higher levels mean a zoom level closer to the Earth's
            # surface with more tiles covering a smaller
            # geographic area and thus more details and vice versa
            # Please note, not all zoom levels are always available
            zoom=14,
            # Use tiles from OpenStreetMap tile server
            source=sourcetiles
        )
        pygmt.makecpt(cmap=cmap, series=lims, background=True)
        fig.grdview(grid=grid, cmap=True, projection=projection, surftype='c', transparency=40)
    else:
        pygmt.makecpt(cmap="gray", series=[-8000, 8000, 1000], continuous=True)
        fig.grdimage(
            grid=topo_data,
            cmap=True,
            region=[minlon, maxlon, minlat, maxlat],
            projection=projection,
            shading=True,
            frame=True
        )
        pygmt.makecpt(cmap=cmap, series=lims, background=True)
        fig.grdimage(grid=grid, cmap=True, projection=projection, frame=True, transparency=40)
    #
    fig.coast(shorelines=True, projection=projection)
    if type(plotvec) != type(None):
        fig.plot(plotvec, projection=projection, region=region)
    fig.colorbar(frame='a10+l"{}"'.format(label))
    # fig.show()
    if interactive:
        return fig, region, projection, xshift, yshift
    else:
        return fig


def get_region(cube):
    '''extracts pyGMT-compliant region from given xr dataset/dataarray'''
    if 'longitude' in cube:
        x='longitude'
        y='latitude'
    elif 'lon' in cube:
        x='lon'
        y='lat'
    elif 'x' in cube:
        x='x'
        y='y'
    else:
        print('no proper dimension')
        return False
    a=cube
    llcrnrlon=float(a[x].min().values) # lower left corner longitude 
    llcrnrlat=float(a[y].min().values) # lower left corner latitude
    urcrnrlon=float(a[x].max().values) # upper right corner longitude
    urcrnrlat=float(a[y].max().values) # upper right corner latitude
    region=[llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat]
    return region


def vis_tif(tifile, stdscale = 1, to_amp_db = False):
    ''' to show a tif file
    Args:
        stdscale (float):  how many stddev are used to form min/max limits for colourscale
        to_amp_db (bool):  will convert intensity to amplitude in dB - useful for LiCSAR epoch mli tifs
    '''
    import rioxarray
    arr = rioxarray.open_rasterio(tifile)
    arr=arr.where(arr !=0)
    if to_amp_db:
        arr.values = np.sqrt(arr.values)
        arr.values = np.log10(arr.values)
    tmean = float(arr.mean())
    tstd = float(arr.std())
    arr.plot(vmin=tmean-stdscale*tstd, vmax=tmean+stdscale*tstd)
    plt.show()


def plot_ts_simple(cube, lon, lat, dvarname = 'cum'):
    ''' returns a simple pygmt figure of time series of xr datacube at given lon, lat coord '''
    time = cube.time.values
    mindate = pd.Timestamp(cube.time.min().values).to_pydatetime().date()
    maxdate = pd.Timestamp(cube.time.max().values).to_pydatetime().date()
    cum_point = cube[dvarname].sel(lon=lon, lat=lat, method='nearest')
    #if dvarname == 'amplitude':
    #    cum_point = 20 * np.log10(cum_point.values)
    #
    if dvarname.find('coh') > -1:
        miny, maxy = 0, 1
    else:
        #miny = float(cube[dvarname].min())
        #maxy = float(cube[dvarname].max())
        miny = float(cum_point.min())
        maxy = float(cum_point.max())
        # stddev=float(cube[dvarname].std())
        # mean=float(cube[dvarname].mean())
        # miny = mean-2*stddev
        # maxy = mean+2*stddev
    #
    fig = pygmt.Figure()
    #
    fig.plot(
        projection="X9c/4.5c",
        region=[mindate, maxdate, miny, maxy],
        # frame=["+t time series (%.3f, %.3f)" % (lon, lat),\
        frame=["xa1Yfg1Y", "yafg+l" + label],
        # "xa1Yfg1Y", "ya10f5+ldisplacement [mm]"],
        x=time,
        y=cum_point,
        style="c0.1c",
        fill="black"
    )
    '''
    if dvarname2:
        fig.plot(
            x=time,
            y=cube[dvarname2].sel(lon=lon, lat=lat, method='nearest'),
            style="c0.1c",
            fill="red"
        )
    '''
    fig.basemap(frame=True)
    return fig

