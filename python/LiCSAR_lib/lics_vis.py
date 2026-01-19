import pygmt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
try:
    import framecare as fc
except:
    print('error loading framecare library - some LiCSAR-connected functionality will not work')

try:
    import contextily as ctx
except:
    print('contextily is not installed')

def plot3(A,B,C, unit='rad', cmap='RdBu', minmax=None):
    '''inputs are three xr.dataarrays to plot

    minmax can be e.g. [-4,4] '''
    origfigsize = plt.rcParams['figure.figsize']
    plt.rcParams["figure.figsize"] = [18,4]
    plt.subplot(1,3,1)
    #AA.origpha
    if minmax:
        A.rename(unit).plot(cmap=cmap, vmin=minmax[0], vmax=minmax[1])
    else:
        A.rename(unit).plot()
    plt.subplot(1,3,2)
    #AA.unwlow
    if minmax:
        B.rename(unit).plot(cmap=cmap, vmin=minmax[0], vmax=minmax[1])
    else:
        B.rename(unit).plot()
    plt.subplot(1,3,3)
    #AA.toremove
    if minmax:
        C.rename(unit).plot(cmap=cmap, vmin=minmax[0], vmax=minmax[1])
    else:
        C.rename(unit).plot()
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



def volcano_clip_plot(volcid, bevel = 0.1, outpng=None):
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
    pygmt.config(FORMAT_GEO_MAP="ddd.xx")
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
        zoom=11,
        # Use tiles from OpenStreetMap tile server
        source=sourcetiles
    )
    # fig.basemap(region=region, projection=projection, frame=["af", '+t"{0}"'.format(title)])
    fig.coast(region=region, shorelines=True, frame=["af"] )  # , water="lightblue")
    fig.plot(volclip.geom, pen='4p,red')
    fig.plot(volcrecord.geom, color='red')
    if outpng:
        fig.savefig(outpng, dpi=150)
    # fig.plot(gpd_overlaps)
    return fig #fig.show()


def pygmt_plot_interactive(cube, title, label='deformation rate [mm/year]', lims=[-15, 50],
                           cmap="polar", photobg=False, plotvec=None, yrange = None):
    ''' This will start a simple interactive viewer in jupyter ntb.
    Note you need ipympl installed and the matplotlib widget must be set here.

    The author is Rochelle Pun, an outstanding student that was selected in highly competitive
    COMET Summer MSc Internship 2024 - and she did magnificent job interactivizing pygmt.
    Her original tool resides here: https://github.com/chelle0425/IntPyGMT

    We only modified this for CIW 2024 tutorial purposes on LiCSBAS

    yrange can be e.g. [-50,50]
    '''
    # print("%matplotlib widget")
    from PIL import Image
    import matplotlib.image as mpimg
    from pygmt.clib import Session
    from pygmt.helpers import GMTTempFile
    #
    # first of all generate the left plot:
    tempng = '/tmp/pygmt_pi.png'
    grid = cube['vel']
    fig, region, projection, xshift, yshift = pygmt_plot(grid, title, label, lims,
                                                         cmap, photobg, plotvec, interactive=True)
    try:
        fig.plot(x=float(cube.ref_lon), y=float(cube.ref_lat), pen="4p,magenta4")
    except:
        print('error adding ref point')
    #
    fig.savefig(tempng)
    #
    png_path = tempng
    spam_path = tempng + "ts.png"
    #
    #time_ds = cube["time"].values
    ddpi = 150
    #
    if isinstance(region, str):
        region_str = region
    elif isinstance(region, list):
        region_str = f"{region[0]}/{region[1]}/{region[2]}/{region[3]}"
    #
    def unpack_xyshift(str):
        if str[-1].isalpha():
            shift_value = str[:-1]
            shift_unit = str[-1]
        else:
            raise Exception("invalid xshift or yshift input")
        return float(shift_value), shift_unit
    #
    xshift_value, xshift_unit = unpack_xyshift(xshift)
    yshift_value, yshift_unit = unpack_xyshift(yshift)
    #
    # we are working in cm
    if xshift_unit == "c":
        xshift_value = xshift_value
    elif xshift_unit == "i":
        xshift_value = xshift_value * 2.54
    elif xshift_unit == "p":
        xshift_value = (xshift_value * 72) * 2.54
    else:
        raise Exception("invalid xshift input (must be either c, i or p)")
    #
    # similarly for y
    if yshift_unit == "c":
        yshift_value = yshift_value
    elif yshift_unit == "i":
        yshift_value = yshift_value * 2.54
    elif yshift_unit == "p":
        yshift_value = (yshift_value * 72) * 2.54
    else:
        raise Exception("invalid yshift input (must be either c, i or p)")
    #
    #    # determine image dimension
    img = Image.open(png_path)
    width, height = img.size  # canvas (width,height) tuple in pixels
    DPI_horz, DPI_vert = img.info.get('dpi')
    #
    #time = pd.to_datetime(time_ds)
    #
    #cum = cube["cum"]
    # ymin = float(cube.cum.mean() - 3*cube.cum.std())
    # ymax = float(cube.cum.mean() + 3*cube.cum.std())
    if yrange:
        ymin,ymax=yrange
    else:
        ymin = float(cube.cum.min())
        ymax = float(cube.cum.max())
    #
    def pos_to_lonlat(x, y):
        # xyshift input in cm
        x = (x / DPI_horz) * 2.54  # convert pixel to cm
        x = x - xshift_value
        y = (y / DPI_vert) * 2.54  # cm
        height_cm = (height / DPI_vert) * 2.54
        y = height_cm - y - yshift_value
        x = [x]  # must be list with one value or np array
        y = [y]
        ### lon lat conversion using mapproject ###
        with Session() as ses:
            with ses.virtualfile_from_vectors(x, y) as fin:
                args = [f'{fin}', f'-R{region_str}', f'-J{projection}', '-I', '-S']
                with GMTTempFile() as fout:
                    ses.call_module(module="mapproject", args=' '.join(args) + " ->" + fout.name)
                    out = fout.read().strip()
        lon, lat = [float(i) for i in out.split(' ')]
        return lon, lat
    #
    def onclick(event):
        # pos.append([event.xdata, event.ydata])
        lon, lat = pos_to_lonlat(event.xdata,
                                 event.ydata)  # pos[-1][0], pos[-1][1]) # pos[-1] represents last click (list with x, y)
        # lonlat.append([lon, lat]) # converts x y to lon lat and appends
        # code below shows how gmt_png can be modified for an interactive time series plot
        # alternatively retrive lon lat using coords_from_figure(ax1)
        # ax1.set_title(f'Click {len(pos)}: {lon}, {lat}')
        fig = plot_ts_simple(cube, lon, lat, label='defo towards sat [mm]', dvarname='cum', miny=ymin, maxy=ymax)
        fig.savefig(spam_path, dpi=ddpi)
        img = mpimg.imread(spam_path)
        ax2.imshow(img, origin='upper')
    #
    # fig = plt.figure(figsize=(12, 5))
    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(121)
    img = mpimg.imread(png_path)
    ax1.imshow(img, origin='upper')
    ax1.axis('off')
    # ax1.set_title(frame+" vel.filt (mm/yr)")
    ax2 = fig.add_subplot(122)
    ax2.axis('off')
    plt.tight_layout()
    #
    label = 'defo towards sat [mm]'
    ffig = plot_ts_simple(cube, float(cube.lon.mean()), float(cube.lat.mean()), label=label, dvarname='cum', miny=ymin,
                          maxy=ymax)
    ffig.savefig(spam_path, dpi=ddpi)
    img = mpimg.imread(spam_path)
    ax2.imshow(img, origin='upper')
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    return ax1


def plotdaz(dazes, frame, toshow = 'daz', lim = 4000, ylim = [-200,200]):
    toplotA = dazes[dazes['AB']=='A'].set_index('epoch')[toshow]*14000
    toplotB = dazes[dazes['AB']=='B'].set_index('epoch')[toshow]*14000
    todelA = toplotA[np.abs(toplotA)>=lim]
    todelB = toplotB[np.abs(toplotB)>=lim]
    print('outlying epochs:')
    if not todelA.empty:
        print(todelA.index.values)
    if not todelB.empty:
        print(todelB.index.values)
    toplotA=toplotA[np.abs(toplotA)<lim]
    toplotB=toplotB[np.abs(toplotB)<lim]
    if not toplotA.empty:
        toplotA.plot(title=frame+' ('+toshow+')',  ylim=ylim, ylabel='$u_{az}$ [mm]', marker='o', color='blue', linestyle='')#-.')
    if not toplotB.empty:
        toplotB.plot(title=frame+' ('+toshow+')',  ylim=ylim, ylabel='$u_{az}$ [mm]', marker='o', color='orange', linestyle='')#-.')
    plt.show()


def plotframedaz(frame, toshow = 'cc_range', lim=4000, ylim1=2000):
    import daz_lib_licsar as dl
    ylim = [-1*ylim1, ylim1]
    dazes = dl.get_daz_frame(frame)
    msab = dl.fc.get_frame_master_s1ab(frame)
    mdatetime=dl.fc.get_master(frame, asdatetime=True)
    epochdates=dazes['epoch'].tolist()
    ABs = dl.flag_s1b(epochdates,mdatetime,msab,True)
    dazes['AB'] = ABs
    plotdaz(dazes, frame, toshow = toshow, lim = lim, ylim=ylim)


## from COPILOT:
import math
import re
def parse_width_cm(code: str) -> int:
    """
    Extract numeric width in cm from a code like 'M13c'.
    Returns the integer width in centimeters.
    """
    match = re.search(r'(\d+)', code)
    if match:
        return int(match.group(1))
    else:
        raise ValueError(f"No numeric width found in '{code}'")


def calculate_zoom(width_cm, dpi, lon_span_deg, latitude_deg=0.0, return_resolution=False):
    """
    Calculate appropriate Web Mercator zoom level for contextily.
    Parameters
    ----------
    width_cm : float
        Width of the map image in centimeters.
    dpi : int
        Resolution in dots per inch.
    lon_span_deg : float
        Longitude span of the map in degrees (WGS-84).
    latitude_deg : float, optional
        Central latitude of the map in degrees (default=0.0, equator).
    Returns
    -------
    int
        Recommended zoom level (0–22).
    """
    # Step 1: convert width to pixels
    inches = width_cm / 2.54
    pixels = inches * dpi
    # Step 2: convert longitude span to meters
    meters_per_degree = 111_320 * math.cos(math.radians(latitude_deg))
    span_meters = lon_span_deg * meters_per_degree
    # Step 3: required resolution (meters per pixel)
    required_res = span_meters / pixels
    if return_resolution:
        return round(np.abs(required_res))
    # Step 4: Web Mercator resolution formula
    # resolution(z) = 156543 / 2^z
    zoom = math.log2(156543 / required_res)
    # Clamp to valid zoom levels
    zoom = max(0, min(22, zoom))
    return round(zoom)


def pygmt_plot(grid, title, label='deformation rate [mm/year]', lims=[-25, 10],
               cmap="roma", photobg=False, plotvec=None, interactive = False,
               region = None, projection = "M13c", medfix=False):
    ''' Function to generate (nice) plot of given grid using pyGMT
    
    Args:
        grid (xr.DataArray): input grid to plot (for tif, use e.g. lics_processing.load_tif2xr)
        title (str):  title (note too long title will disbalance the figure)
        label (str):  label below the colour scale
        lims (list):  colour scale limits (if None, it will do min max minus 2std)
        cmap (str):   colour scale map (try 'vik' for E-W)
        photobg (bool): will plot orthophotomap as the background (if False, DEM relief is used)
        plotvec (geopandas etc): will plot vector data to the map, using pyGMT defaults
        region (tuple/None):  either None or set using tuple (minlon, maxlon, minlat, maxlat)
        projection (str):  standard GMT projection string - try e.g."M13c" or 'R13c' for Robinson etc.
        medfix (bool): if True, it will center the values to their median (set 0)

    Returns:
        pygmt.figure.Figure
    '''
    try:
        grid = grid.load()
        grid = grid.where(grid != 0)
        isgrid = True
    except:
        print('error loading the input dataarray to memory')
        if type(plotvec)==type(None):
            return False
        else:
            print('trying to plot the map using provided vector data (assuming geopandas)')
            region = plotvec.bounds.minx.min(),\
                     plotvec.bounds.maxx.max(),\
                     plotvec.bounds.miny.min(),\
                     plotvec.bounds.maxy.max()
            # x1, x2, y1, y2
            lims = True
            isgrid = False

    # try cmap 'vik' for E-W
    #
    # grid = a['U'].where(a.mask < 5) - 10
    # topo_data = '@earth_relief_03s' #3 arc second global relief (SRTM3S)
    topo_data = '@earth_relief_01s'  # 1 arc second global relief (SRTM3S) (30 m)

    if not region:
        minlon, maxlon = float(np.min(grid.lon)), float(np.max(grid.lon))
        minlat, maxlat = float(np.min(grid.lat)), float(np.max(grid.lat))
    else:
        minlon, maxlon, minlat, maxlat = region

    if not lims:
        tmean = float(grid.mean())
        tstd = float(grid.std())
        stdscale = 2
        lims = [tmean - stdscale * tstd, tmean + stdscale * tstd]

    fig = pygmt.Figure()
    pygmt.config(FORMAT_GEO_MAP="ddd.xx") #, MAP_FRAME_TYPE="plain")
    # projection = "M13c" # 'R13c' for Robinson etc.
    region = [minlon, maxlon, minlat, maxlat]
    if medfix:
        grid=grid - grid.sel(lon=slice(minlon,maxlon), lat=slice(maxlat, minlat)).median()

    if interactive:
        xshift = '1.5c'
        yshift = '2.5c'
        fig.shift_origin(xshift=xshift, yshift=yshift)
    fig.basemap(region=region, projection=projection, frame=["af", '+t"{0}"'.format(title)])

    if photobg:
        import contextily as ctx
        sourcetiles = ctx.providers.Esri.WorldImagery
        try:
            width_cm=parse_width_cm(projection)
            zoomlevel=calculate_zoom(width_cm, dpi=150, lon_span_deg=maxlon-minlon, latitude_deg=round((maxlat-minlat)/2, 1))
            print('using zoom level of '+str(zoomlevel))
        except:
            print('error calculating zoom level - using default')
            zoomlevel=9
        fig.tilemap(
            region=region, projection=projection,
            # region=[-157.84, -157.8, 21.255, 21.285],
            # projection="M12c",
            # Set level of details (0-22)
            # Higher levels mean a zoom level closer to the Earth's
            # surface with more tiles covering a smaller
            # geographic area and thus more details and vice versa
            # Please note, not all zoom levels are always available
            zoom=zoomlevel,
            # Use tiles from OpenStreetMap tile server
            source=sourcetiles
        )
        if isgrid:
            pygmt.makecpt(cmap=cmap, series=lims, background=True)
            fig.grdview(grid=grid, cmap=True, projection=projection, surftype='c', transparency=40)
    else:
        try:
            width_cm=parse_width_cm(projection)
            imresol=calculate_zoom(width_cm, dpi=150, lon_span_deg=maxlon-minlon, latitude_deg=round((maxlat-minlat)/2, 1), return_resolution=True)
            # https://docs.generic-mapping-tools.org/6.0/datasets/earth_relief.html
            if imresol<80:
                topo_data = '@earth_relief_01s'  # 1 arc second global relief (SRTM3S) (30 m)
            elif imresol<400:
                topo_data = '@earth_relief_03s'
            elif imresol<800:
                topo_data = '@earth_relief_15s'
            elif imresol<1700:
                topo_data = '@earth_relief_30s'
            else:
                topo_data = '@earth_relief_01m'
            print('using topo_data '+str(topo_data))
        except:
            print('error calculating im resolution - using default topo_data: '+str(topo_data))
        pygmt.makecpt(cmap="gray", series=[-8000, 8000, 1000], continuous=True)
        fig.grdimage(
            grid=topo_data,
            cmap=True,
            region=[minlon, maxlon, minlat, maxlat],
            projection=projection,
            shading=True,
            frame=True
        )
        if isgrid:
            pygmt.makecpt(cmap=cmap, series=lims, background=True)
            fig.grdimage(grid=grid, cmap=True, projection=projection, frame=True, transparency=40)
    #
    fig.coast(shorelines=True, projection=projection)
    if type(plotvec) != type(None):
        fig.plot(plotvec, projection=projection, region=region)
    # for colorbar, use proper step:
    rng = lims[1] - lims[0]
    step=round_to_nice_step(rng/8)
    if str(step).split('.')[-1] == '0':
        step=int(step)
    #fig.colorbar(frame='a10+l"{}"'.format(label))
    fig.colorbar(frame='a{0}+l"{1}"'.format(str(step),label))
    # fig.show()
    if interactive:
        return fig, region, projection, xshift, yshift
    else:
        return fig


def round_to_nice_step(step):
    # Get the order of magnitude
    magnitude = 10 ** np.floor(np.log10(step))
    # Normalize step to [1, 10)
    normalized = step / magnitude
    # Choose a nice base step
    if normalized < 1.5:
        nice = 1
    elif normalized < 3:
        nice = 2
    elif normalized < 7:
        nice = 5
    else:
        nice = 10
    return nice * magnitude


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


def generate_vel_preview(nc, lims = [-100,100], title = 'vel_msk', volcid = None, outpng = None):
    import geopandas as gpd
    import xarray as xr
    import volcdb as v
    from shapely.geometry import Point
    import framecare as fc
    #
    cube = xr.open_dataset(nc)
    refpoint = Point((cube.ref_lon, cube.ref_lat))
    refpoint = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[refpoint])
    #
    print('Plotting velocity layer on DEM')
    grid = cube.vel.where(cube.mask==1)
    if volcid:
        lbstr = v.get_licsbas_clipstring_volcano(volcid) # size as in volc portal
        minlon, maxlon, minlat, maxlat = [ float(s) for s in lbstr.split('/') ]
        grid = grid.sel(lon=slice(minlon,maxlon), lat=slice(minlat,maxlat)).load()
        volctb = v.get_volc_info(volcid)
        #title = str(volctb.name.values[0])
        poly = fc.lonlat_to_poly(minlon, maxlon, minlat, maxlat)
        allvolcsinpoly = v.get_volcanoes_in_polygon(poly, volcs=None)
    #
    gridgmt = pygmt_plot(grid, title, 'mm/year', lims, 'vik', region = (minlon, maxlon, minlat, maxlat))
    gridgmt.plot(data=refpoint, pen="4p,magenta4")
    if volcid:
        gridgmt.plot(data=allvolcsinpoly.geom, style="t6p", pen="0.5p,black", fill="blue")
        gridgmt.plot(data=volctb.geom, style="t8p", pen="0.5p,black", fill="red")
    if outpng:
        gridgmt.savefig(outpng, dpi=120)
    return gridgmt


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


def plot_ts_simple(cube, lon, lat, label = 'test', dvarname = 'cum', miny=None, maxy=None):
    ''' returns a simple pygmt figure of time series of xr datacube at given lon, lat coord '''
    time = cube.time.values
    mindate = pd.Timestamp(cube.time.min().values).to_pydatetime().date()
    maxdate = pd.Timestamp(cube.time.max().values).to_pydatetime().date()
    cum_point = cube[dvarname].sel(lon=lon, lat=lat, method='nearest')
    if type(miny)==type(None):
        if dvarname.find('coh') > -1:
            miny, maxy = 0, 1
        else:
            miny = float(cum_point.min())
            maxy = float(cum_point.max())
    #
    fig = pygmt.Figure()
    pygmt.config(FONT_LABEL="10p")
    #
    with pygmt.config(FONT_TITLE="12p"):
        fig.plot(
            projection="X9c/7c",
            region=[mindate, maxdate, miny, maxy],
            # frame=["+t time series (%.3f, %.3f)" % (lon, lat),\
            frame=["xa1Yfg1Y", "yafg+l" + label, "+t%.3f, %.3f" % (lon, lat)],
            # "xa1Yfg1Y", "ya10f5+ldisplacement [mm]"],
            x=time,
            y=cum_point,
            style="c0.1c",
            fill="black"
        )
    fig.basemap(frame=True)
    return fig


def pygmt_plot_bursts(burstlist = None, framelist = None, title = '', background = True,
                      projection = "M8c", bevel = 0.05, label_nofiles = False, label_bids = True):
    sourcetiles = ctx.providers.Esri.WorldImagery
    if framelist:
        framesgpd = fc.get_frames_gpd(framelist)
    #
    if burstlist:
        bidsgpd = fc.bursts2geopandas(burstlist)
        # framelist = ['037D_04691_021207']
        region = bidsgpd.bounds.minx.min() - bevel, \
                 bidsgpd.bounds.maxx.max() + bevel, \
                 bidsgpd.bounds.miny.min() - bevel, \
                 bidsgpd.bounds.maxy.max() + bevel
    elif framelist:
        region = framesgpd.bounds.minx.min() - bevel, \
                 framesgpd.bounds.maxx.max() + bevel, \
                 framesgpd.bounds.miny.min() - bevel, \
                 framesgpd.bounds.maxy.max() + bevel
    else:
        print('nothing to plot - provide burst or framelist')
        return False
    #
    fig = pygmt.Figure()
    # pygmt.config(FORMAT_GEO_MAP="ddd.x", MAP_FRAME_TYPE="plain")
    pygmt.config(FONT_TITLE="12p")
    pygmt.config(FONT_LABEL="8p")
    if title:
        fig.basemap(region=region, projection=projection, frame=["af", '+t"{0}"'.format(title)])
    else:
        fig.basemap(region=region, projection=projection, frame=["af"])
    #
    if background:
        fig.tilemap(
            region=region,
            zoom=8,
            # Use tiles from OpenStreetMap tile server
            source=sourcetiles
        )
    #
    fig.coast(region=region, shorelines=True, frame=True)
    if burstlist:
        fig.plot(bidsgpd.geometry, pen='1p,red')
    if framelist:
        # framesgpd = fc.get_frames_gpd(framelist)
        fig.plot(framesgpd.geometry, pen='1p,blue')
        fig.text(text=framelist, x=framesgpd.geometry.centroid.x.values,
                 y=framesgpd.geometry.centroid.y.values,
                 font="8p,Helvetica,blue")
    #
    if burstlist:
        if label_nofiles:
            iwbs = fc.bursts_group_to_iws(burstlist)
            iwl = 0
            i = 0
            for iw in iwbs:
                iwlt = len(iw)
                if iwlt > iwl:
                    iwl = iwlt
                    seliw = i
                i += 1
            iwbs = iwbs[seliw]
            bidsgpd = fc.bursts2geopandas(iwbs)
            blens = []
            print('getting no of files per burst, as imported to LiCSInfo')
            for b in iwbs:
                blens.append(len(fc.get_files_from_burst(b)))
            bidsgpd['nofiles'] = blens
            fig.text(text=bidsgpd.nofiles.values, x=bidsgpd.geometry.centroid.x.values,
                     y=bidsgpd.geometry.centroid.y.values,
                     font="8p,Helvetica")
        elif label_bids:
            fig.text(text=bidsgpd.burstID.values, x=bidsgpd.geometry.centroid.x.values,
                     y=bidsgpd.geometry.centroid.y.values,
                     font="5p,Helvetica")  # , justify="RT", offset='J/30p')
    #
    return fig