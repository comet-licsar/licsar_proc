#!/usr/bin/env python
"""
Tool which levereges the LiCSAR database to convert a polygon file to lists of 
files stored on tape (NLA) and on disk.

usage:

    polygon2filelist -p polygonfile -d ondisklist -t ontapelist
"""

################################################################################
# imports
################################################################################

import LiCSquery as lq
import re
import os.path as path
import sys
import copy
from getopt import getopt
from fastkml import kml
from fastkml import styles
from shapely.geometry import Polygon
################################################################################
# Constants
################################################################################
coordRE="(-?\d+\.\d+)\s?[,-]?\s?(-?\d+\.\d+)" # Coord regualar expresion 


################################################################################
# get frames which intersect (lazy) the polygon
################################################################################
def get_frames_in_polygon(lon1,lat1,lon2,lat2):
    sql_q = "select polyid,polyid_name, active,"\
            'corner1_lon, corner1_lat, '\
            'corner2_lon, corner2_lat, '\
            'corner3_lon, corner3_lat, '\
            'corner4_lon, corner4_lat, '\
            'corner5_lon, corner5_lat, '\
            'corner6_lon, corner6_lat, '\
            'corner7_lon, corner7_lat, '\
            'corner8_lon, corner8_lat, '\
            'corner9_lon, corner9_lat, '\
            'corner10_lon, corner10_lat, '\
            'corner11_lon, corner11_lat, '\
            'corner12_lon, corner12_lat '\
            "from polygs " \
            "where "\
            "greatest( "\
            "corner1_lon,corner2_lon,corner3_lon,corner4_lon,"\
            "corner5_lon,corner6_lon,corner7_lon,corner8_lon,"\
            "corner9_lon,corner10_lon,corner11_lon,corner12_lon"\
            ") "\
            ">= {0} and "\
            "least( "\
            "corner1_lon,corner2_lon,corner3_lon,corner4_lon,"\
            "corner5_lon,corner6_lon,corner7_lon,corner8_lon,"\
            "corner9_lon,corner10_lon,corner11_lon,corner12_lon"\
            ") "\
            "<= {1} ".format(lon1,lon2)
    sql_q += "and "\
            "greatest( "\
            "corner1_lat,corner2_lat,corner3_lat,corner4_lat,"\
            "corner5_lat,corner6_lat,corner7_lat,corner8_lat,"\
            "corner9_lat,corner10_lat,corner11_lat,corner12_lat"\
            ") "\
            ">= {0} and "\
            "least( "\
            "corner1_lat,corner2_lat,corner3_lat,corner4_lat,"\
            "corner5_lat,corner6_lat,corner7_lat,corner8_lat,"\
            "corner9_lat,corner10_lat,corner11_lat,corner12_lat"\
            ") "\
            "<= {1};".format(lat1,lat2)
    return lq.do_query(sql_q)

def get_bursts_in_polygon(lon1,lat1,lon2,lat2):
    sql_q = "select bid ,bid_tanx, "\
            'corner1_lon, corner1_lat, '\
            'corner2_lon, corner2_lat, '\
            'corner3_lon, corner3_lat, '\
            'corner4_lon, corner4_lat '\
            "from bursts " \
            "where "\
            "greatest( "\
            "corner1_lon,corner2_lon,corner3_lon,corner4_lon"\
            ") "\
            ">= {0} and "\
            "least( "\
            "corner1_lon,corner2_lon,corner3_lon,corner4_lon"\
            ") "\
            "<= {1} ".format(lon1,lon2)
    sql_q += "and "\
            "greatest( "\
            "corner1_lat,corner2_lat,corner3_lat,corner4_lat"\
            ") "\
            ">= {0} and "\
            "least( "\
            "corner1_lat,corner2_lat,corner3_lat,corner4_lat"\
            ") "\
            "<= {1};".format(lat1,lat2)
    return lq.do_query(sql_q)

################################################################################
# Create KML styles
################################################################################
def createKMLStyles(colours,parent):
    lnStyle = styles.LineStyle(color='FFFFFFFF')
    polStyle = styles.PolyStyle(color='FFFFFFFF')
    styl = styles.Style(id='default')
    styl.append_style(lnStyle)
    styl.append_style(polStyle)
    parent.append_style(styl)

    for cName,cVal in colours.iteritems():
        styl = styles.Style(id=cName)
        styl.append_style(lnStyle)
        polStyle = styles.PolyStyle(color=cVal)
        styl.append_style(polStyle)
        parent.append_style(styl)
    
################################################################################
# Create KML for polygons
################################################################################
def createFrameKML(polygonName,polygon,frames,colours,parent):
        
    folder = kml.Folder(name=polygonName+' Frames',description='LiCSAR Frames within '+polygonName)

    
    stylNames = colours.keys()
    styleIter = 0

    plceMrk = kml.Placemark(name=polygonName,description='Master Polygon',styleUrl='default')
    plceMrk.geometry = polygon

    folder.append(plceMrk)

    for frame in frames:
        frmPoly = Polygon([[frame[i],frame[i+1]] for i in range(3,27,2)])
        if frmPoly.intersects(polygon):
            plceMrk = kml.Placemark(name=frame[1],styleUrl=stylNames[styleIter])
            styleIter = (styleIter+1)%len(stylNames)
            plceMrk.geometry = Polygon(frmPoly)
            folder.append(plceMrk)

    parent.append(folder)
    return parent

def createBurstKML(polygonName,polygon,bursts,colours,parent):
        
    folder = kml.Folder(name=polygonName+' Bursts',description='LiCSAR Bursts within '+polygonName)

    
    stylNames = colours.keys()
    styleIter = 0

    plceMrk = kml.Placemark(name=polygonName,description='Master Polygon',styleUrl='default')
    plceMrk.geometry = polygon

    folder.append(plceMrk)

    for burst in bursts:
        brstPoly = Polygon([[burst[i],burst[i+1]] for i in range(2,10,2)])
        if brstPoly.intersects(polygon):
            plceMrk = kml.Placemark(name=burst[1],styleUrl=stylNames[styleIter])
            styleIter = (styleIter+1)%len(stylNames)
            plceMrk.geometry = Polygon(brstPoly)
            folder.append(plceMrk)

    parent.append(folder)
    return parent

################################################################################
# Create KML from Source KML
################################################################################
def createFrameKMLFromSourceKML(sourceKMLFile,frameKMLFile):
    polygs = get_polys_from_kml(sourceKMLFile)
    baseKML = kml.KML()
    doc = kml.Document(name='Frames from polygons',description='LiCSAR frames within polygons')

    colours = {'red':'FF0000FF',
            'green':'00FF00FF',
            'blue':'0000FFFF',
            'purple':'CC00FFFF'}

    createKMLStyles(colours,doc)
    
    for poly in polygs:
        polyName = poly[0]
        polyBnds = poly[1].bounds
        frms = get_frames_in_polygon(*polyBnds)
        createFrameKML(polyName,poly[1],frms,colours,doc)

    baseKML.append(doc)

    with open(frameKMLFile,'w') as f:
        f.write(baseKML.to_string())

def createBurstKMLFromSourceKML(sourceKMLFile,burstKMLFile):
    polygs = get_polys_from_kml(sourceKMLFile)
    baseKML = kml.KML()
    doc = kml.Document(name='Frames from polygons',description='LiCSAR bursts within polygons')

    colours = {'red':'FF0000FF',
            'green':'00FF00FF',
            'blue':'0000FFFF',
            'purple':'CC00FFFF'}

    createKMLStyles(colours,doc)
    
    for poly in polygs:
        polyName = poly[0]
        polyBnds = poly[1].bounds
        brsts = get_bursts_in_polygon(*polyBnds)
        createBurstKML(polyName,poly[1],brsts,colours,doc)

    baseKML.append(doc)

    with open(burstKMLFile,'w') as f:
        f.write(baseKML.to_string())
################################################################################
# Get polygons (shapely) from kml
################################################################################
def get_polys_from_kml(kmlfile):

    def parse_kml_tree(k):
        if type(k) == kml.Placemark:
            return [[k.name,k.geometry]]
        else:
            polys = list()
            for ftr in k.features():
                polys += parse_kml_tree(ftr)
            return polys

    with open(kmlfile,'r') as f:
        kmlDoc = f.read()

    baseKml = kml.KML()
    baseKml.from_string(kmlDoc)

    return parse_kml_tree(baseKml)


################################################################################
# Get polygon from file
################################################################################
def getPolygonFromFile(polyFile):
    
    # Setup polygon set
    polygon = list()

    #iterate over polygon file
    with open(polyFile) as pf:
        for line in pf:
            mtch = re.search(coordRE,line)
            coord = [float(numStr) for numStr in mtch.groups()]
            polygon.append(coord)

    return polygon

################################################################################
# Get polygon long/lat min max
################################################################################
def getPolygonMinMax(polygon):

    polyLon = [coord[0] for coord in polygon]
    polyLat = [coord[1] for coord in polygon]

    return min(polyLon),max(polyLon),min(polyLat),max(polyLat)

################################################################################
# Convert polygon limits to file list
################################################################################
def getFilesFromPolygonMinMax(polyMinMax,track):

    files = set()
    brstTracks = set()

    brstIDs = [brst[0] for brst in lq.get_bursts_in_polygon(*polyMinMax)]

    for brstID in brstIDs:
        curTrack = int(re.search('(\d+)_',brstID).groups()[0])
        if not track:
            brstTracks.add(curTrack)
        brstFiles = [bf[1] for bf in lq.get_files_from_burst(brstID)]
        if curTrack==track:
            files |= set(brstFiles)

    return list(files),list(brstTracks)
