import geopandas as gpd
import folium
import volcdb as v
import requests
from bs4 import BeautifulSoup
import os

print('see /gws/ssde/j25a/nceo_geohazards/vol1/public/shared/temp/earmla/volcano_map')
url = "https://comet-volcanodb.org/testing/licsbas"

# first, extract available volcanoes
def get_web_volcanoes(url, return_sublinks = True):
    html = requests.get(url).text
    soup = BeautifulSoup(html, "html.parser")
    # Extract all <li><a>...</a></li> text
    if return_sublinks:
        volcano_sublinks = [a['href'].strip() for a in soup.select("ul li a")]
        return volcano_sublinks
    else:
        volcano_names = [a.text.strip() for a in soup.select("ul li a")]
        return volcano_names

volclinks = get_web_volcanoes()
volcpd = v.get_volc_info()

# get selection from the pd:
volcnames= []
volcgeom = []
for vl in volclinks:
    vlpd = volcpd[volcpd['vportal_name']==vl.lower().split('.')[0]]
    skip = False
    if vlpd.empty:
        # try find:
        vname = vl.split('.')[0].split('_')[0]
        vlpd = v.find_volcano_by_name(vname)
        if len(vlpd)>1:
            #print('more finds from '+vname)
            vlpd = v.find_volcano_by_name(vl.split('.')[0].split('_')[1])
        if len(vlpd)>1:
            print('more finds from '+vname)
            print('skipping')
            skip = True
        elif vlpd.empty:
            print('not found - skipping')
            skip = True
    if skip:
        volcnames.append(None)
        volcgeom.append(None)
    else:
        volcnames.append(vlpd['name'].values[0])
        volcgeom.append(vlpd['geom'].values[0])


import geopandas as gpd
#from shapely import wkt

# Convert WKT strings → Shapely geometry objects
#geoms = [wkt.loads(g) for g in volcgeom]

# Build GeoDataFrame
gdf = gpd.GeoDataFrame(
    {
        "name": volcnames,
        "link": [os.path.join(url, a) for a in volclinks]
    },
    geometry=volcgeom,
    crs="EPSG:4326"   # assume WGS84; change if needed
)

# removing not found:
gdf = gdf.dropna()

# Your GeoDataFrame: gdf (with columns: vportal_name, geometry)
# Make sure it's in WGS84
# gdf = gdf.to_crs(4326)

# Create folium map centered roughly at your data
m = folium.Map(location=[0, 20], zoom_start=2)

for _, row in gdf.iterrows():
    lat = row.geometry.y
    lon = row.geometry.x
    name = row["name"]
    licsbasurl = row['link']  # f"https://web.page/{name}"
    licsalerturl = licsbasurl.replace('/licsbas/', '/licsalert/')
    popup_html = f"{name}:<br /><a href='{licsbasurl}' target='_blank'>LiCSBAS</a><br /><a href='{licsalerturl}' target='_blank'>LiCSAlert</a>"

    folium.Marker(
        location=[lat, lon],
        popup=popup_html,
        tooltip=name,
        # icon=folium.Icon(color='red', icon='glyphicon-fire')
        icon=folium.DivIcon(
            html=f"""<div style="font-size: 20px; color: red;">&#9650;</div>""")
    ).add_to(m)

title_html = """
     <h3 align="center" style="font-size:20px"><b>DEEPVOLC Pilot Volcano Deformation Portal</b></h3>
     """

m.get_root().html.add_child(folium.Element(title_html))

# Save to HTML file
m.save("volcano_map.html")


