#!/usr/bin/env python

import geopandas as gpd
import folium
import volcdb as v
import requests
from bs4 import BeautifulSoup
import os
from branca.element import Element

checknisar=True
print('see /gws/ssde/j25a/nceo_geohazards/vol1/public/shared/temp/earmla/volcano_map')
url = "https://comet-volcanodb.org/testing/licsbas"
outhtml = '/gws/ssde/j25a/nceo_geohazards/vol1/public/shared/temp/earmla/volcano_map/volcano_map.html'

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

if checknisar:
    import nisardata as nd


volclinks = get_web_volcanoes(url)
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
            if vl.split('.')[0].split('_')[1] == 'Nevada':
                print('assuming Sierra Nevada - volcano ID 355123')
                vlpd = vlpd[vlpd.volc_id == 355123]
            else:
                print('skipping')
                skip = True
        elif vlpd.empty:
            print('not found - skipping')
            skip = True
    if skip:
        volcnames.append(None)
        volcgeom.append(None)
    else:
        if checknisar:
            volcanoid = vlpd.volc_id.values[0]
            nisars = nd.get_nisar_data_for_volcano(volcanoid)
            nislen = len(nisars)
            if nislen > 0:
                print('GREAT NEWS - THERE ARE '+str(nislen)+' NISAR GSLCs for volcid '+str(volcanoid)+' - contact Milan')
        volcnames.append(vlpd['name'].values[0])
        volcgeom.append(vlpd['geom'].values[0])


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

# footer:


footer_html = """
<style>
/** footer: **/

.footer {
  color: #ffffff;
  background-color: #272727;
  min-height: 2em;
}

.footer a {
  color: #ffffff;
  text-decoration: none;
  outline: none;
  padding-left: 5px;
  padding-right: 5px;
}

.footer_text {
  display: flex;
  flex: 1 1 auto;
  flex-flow: row wrap;
  justify-content: center;
  text-align: center;
  padding-top: 0.1em;
  padding-bottom: 0em;
  padding-left: 1.5em;
  padding-right: 1.5em;
  color: #dddddd;
  background-color: #272727;
}

.footer_text a {
  color: #ffffff;
  outline: none;
}

.footer_logos {
  display: flex;
  flex: 1 1 auto;
  flex-flow: row wrap;
  justify-content: center;
  background-color: #272727;
}

.footer_item {
  display: inline-block;
  padding-bottom: 0.1em;
  background-color: #272727;
}

.footer_logo {
  -webkit-filter: grayscale(100%);
  filter: grayscale(100%);
  max-height: 4em;
  max-width: 270px;
  margin: 0.5em;
  background-color: #272727;
}

/* Make footer stick to bottom */
#footer {
    position: fixed;
    bottom: 0;
    left: 0;
    width: 100%;
    z-index: 9999;
}

/* Prevent footer covering parts of the map */
.leaflet-container {
    padding-bottom: 90px;
}

</style>

<div id="footer" class="footer">
    <div id="footer_text" class="footer_text">
        <p>funded by DeepVolc, COMET, and ERC</p>
    </div>

    <div id="footer_logos" class="footer_logos">
        <div id="footer_deepvolc_logo" class="footer_item">
            <a href="https://environment.leeds.ac.uk/dir-record/research-projects/1801/forecasting-volcanic-activity-using-deep-learning-deepvolc/"
               target="_blank">
              <img class="footer_logo" src="https://comet-volcanodb.org/testing/licsbas/img/footer/deepvolc_logo.png">
            </a>
        </div>

        <div id="footer_comet_logo" class="footer_item">
            <a href="https://comet.nerc.ac.uk/"
               target="_blank">
              <img class="footer_logo" src="https://comet-volcanodb.org/testing/licsbas/img/footer/comet_logo.png">
            </a>
        </div>

        <div id="footer_erc_logo" class="footer_item">
            <a href="https://erc.europa.eu/"
               target="_blank">
              <img class="footer_logo" src="https://comet-volcanodb.org/testing/licsbas/img/footer/erc_logo.png">
            </a>
        </div>
    </div>
</div>
"""

m.get_root().html.add_child(Element(footer_html))

# Save to HTML file
m.save(outhtml)


