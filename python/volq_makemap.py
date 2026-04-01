import geopandas as gpd
import folium

# Your GeoDataFrame: gdf (with columns: vportal_name, geometry)
# Make sure it's in WGS84
gdf = gdf.to_crs(4326)

# Create folium map centered roughly at your data
m = folium.Map(location=[0, 20], zoom_start=2)

for _, row in gdf.iterrows():
    lat = row.geom.y
    lon = row.geom.x
    name = row["vportal_name"]
    url = f"https://web.page/{name}"

    popup_html = f"<a href='{url}' target='_blank'>{name}</a>"

    folium.Marker(
        location=[lat, lon],
        popup=popup_html,
        tooltip=name
    ).add_to(m)

# Save to HTML file
m.save("volcano_map.html")
``


#### or
import folium
import geopandas as gpd
import json

gdf = gdf.to_crs(4326)

# Add a URL field
gdf["url"] = "https://web.page/" + gdf["vportal_name"]

# Convert to geojson
geojson = gdf.to_json()

m = folium.Map(location=[0,20], zoom_start=2)

folium.GeoJson(
    geojson,
    name="volcanoes",
    popup=folium.GeoJsonPopup(fields=["vportal_name", "url"]),
    tooltip=folium.GeoJsonTooltip(fields=["vportal_name"])
).add_to(m)

m.save("volcano_map.html")
