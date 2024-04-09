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


def volcano_clip_plot(volcid, bevel = 0.1):
    '''plots volcano given by its volcano id - see volcdb'''
    #volcid = 243040
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


def pygmt_plot(grid, title, label='deformation rate [mm/year]', lims=[-25, 10],
               cmap="roma", photobg=False, plotvec=None):
    # try cmap 'vik' for E-W
    #
    # grid = a['U'].where(a.mask < 5) - 10
    # topo_data = '@earth_relief_03s' #3 arc second global relief (SRTM3S)
    topo_data = '@earth_relief_01s'  # 3 arc second global relief (SRTM3S)

    minlon, maxlon = float(np.min(grid.lon)), float(np.max(grid.lon))
    minlat, maxlat = float(np.min(grid.lat)), float(np.max(grid.lat))

    fig = pygmt.Figure()
    projection = "M13c" # 'R13c' for Robinson etc.
    region = [minlon, maxlon, minlat, maxlat]
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
    return fig


def vis_tif(tifile, stdscale = 1, todb = False):
    import rioxarray
    arr = rioxarray.open_rasterio(tifile)
    arr=arr.where(arr !=0)
    if todb:
        arr.values = np.log10(arr.values)
    tmean = float(arr.mean())
    tstd = float(arr.std())
    arr.plot(vmin=tmean-stdscale*tstd, vmax=tmean+stdscale*tstd)
    plt.show()
