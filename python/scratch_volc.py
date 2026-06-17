import volcdb as v
from shapely.geometry import Polygon, box
import geopandas as gpd
import lics_vis as lv
import os

def enlarge_rectangle(polygon, distance, cstr='RLUD'):
    minx, miny, maxx, maxy = polygon.bounds
    new_minx = minx #- distance
    new_miny = miny #- distance
    new_maxx = maxx # + distance
    new_maxy = maxy
    if 'L' in cstr:
        new_minx = minx - distance
    if 'D' in cstr:
        new_miny = miny - distance
    if 'R' in cstr:
        new_maxx = maxx + distance
    if 'U' in cstr:
        new_maxy = maxy + distance
    return Polygon([(new_minx, new_miny), (new_minx, new_maxy), (new_maxx, new_maxy), (new_maxx, new_miny)])


# merging with other volcano
def merge_volclips_of_volcanoes(volcid, vtomerge):
    vidclip = v.get_volclip_vids(volcid)[0]
    mergeclipid = v.get_volclip_vids(vtomerge)[0]
    vidclip_geom = v.get_volclips_gpd(vidclip).geom.values[0]
    vtomerge_geom = v.get_volclips_gpd(mergeclipid).geom.values[0]
    return merge_gpds(vidclip_geom, vtomerge_geom)


def merge_gpds(vidclip_geom, vtomerge_geom):
    minx, miny, maxx, maxy = vtomerge_geom.union(vidclip_geom).bounds
    newpoly = box(minx, miny, maxx, maxy)
    newpoly = gpd.GeoDataFrame({'geom': [newpoly]})  # , crs=largerpr.crs)
    newpoly = newpoly.set_geometry('geom')
    return newpoly


############ update volcano clip - Pedro's boxes
# 1. identify the volcano and set the preliminary box
# 2. plot the situation
# 3. when all fine, run the re-volclip procedure

### input data
# volcano ID
vid = 241070
# vname = 'Taupo'
# Pedro's clip
largeclip = True #False
new_minx, new_maxx, new_miny, new_maxy = 119.9528, 120.5624, 15.3599, 14.954
#    175.498, 176.315, -39.055, -38.536
vtomerge = None

#### 1, and 2,
# get clip poly
vclipid = v.get_volclip_vids(vid)[0]
origpoly=v.get_volclips_gpd(vclipid).geom.values[0]
# transform to standard 50x50 km
distorig = (origpoly.bounds[2] - origpoly.bounds[0])*111.111
bufferdistance = (50 - distorig)/2/111.111
newpoly=enlarge_rectangle(origpoly, bufferdistance)
#newpoly = enlarge_rectangle(newpoly, -4 / 111.111, cstr='D')
#newpoly = enlarge_rectangle(newpoly, 4 / 111.111, cstr='L')
gg=gpd.GeoDataFrame({'geom': [newpoly]})
gg=gg.set_geometry('geom') #.values
# transform Pedro's clip
if largeclip:
    largepoly = Polygon([(new_minx, new_miny), (new_minx, new_maxy), (new_maxx, new_maxy), (new_maxx, new_miny)])
    #largepoly = enlarge_rectangle(largepoly, -20 / 111.111, cstr='LR')
    largepoly = enlarge_rectangle(largepoly, -7 / 111.111, cstr='L')
    #largepoly = enlarge_rectangle(largepoly, -8 / 111.111, cstr='UD')
    largepoly = gpd.GeoDataFrame({'geom': [largepoly]}) #, crs=largerpr.crs)
    largepoly=largepoly.set_geometry('geom')

# or if merging:
if vtomerge:
    newpoly = merge_volclips_of_volcanoes(vid, vtomerge)
    #gg = enlarge_rectangle(newpoly.geom.values[0], 8/111.111, cstr='D') # UDLR
    gg = enlarge_rectangle(newpoly.geom.values[0], 5/111.111, cstr='DR')
    gg = enlarge_rectangle(gg, 2/111.111, cstr='U')
    gg=gpd.GeoDataFrame({'geom': [gg]})
    gg=gg.set_geometry('geom') #.values

# plot the overview
fig=lv.volcano_clip_plot_with_frames(vid)
fig.plot(gg.geom, pen='0.5p,orange')
if largeclip:
    fig.plot(largepoly.geom, pen='1p,red')

fig.savefig('delme.png', dpi=150); os.system('display delme.png')



#### 3,
# if ok, change things:
newpoly = gg   # if 50x50 km (or merging)
# newpoly = largepoly  # if Pedro's clip
primvolc = vid
lon1, lat1, lon2, lat2 = newpoly.geom[0].bounds
new_cliparea_geo = str(lon1) + '/' + str(lon2) + '/' + str(lat1) + '/' + str(lat2)
print('old clip ID: '+str(vclipid))

# deleting volclip from db, links from procdir, and moving to backup dir temporarily
vclipid = v.get_volclip_vids(vid)[0]
v.delete_volclip(vclipid)
import glob
oldfrs = glob.glob(os.path.join(os.environ['LiCSAR_volc'], str(vclipid))+'/*/corner*')
for old in oldfrs:
    fr=old.split('.')[-1]
    tr=str(int(fr[:3]))
    todel=os.path.join(os.environ['LiCSAR_procdir'], tr, fr, 'subsets', str(vclipid))
    os.system('rm '+todel)


cmd = "mv $LiCSAR_volc/{0} $LiCSAR_volc/../volc_subsets.backup/{1}".format(str(vclipid), str(vid))
os.system(cmd)
# os.system("echo volcid="+str(vid)+" > $LiCSAR_volc/../volc_subsets.backup/"+str(vclipid)+"/volcid")


newclip = v.create_volclip_for_volcano(vid, new_cliparea_geo)
print('new clip ID: '+str(newclip))

# changing for volcanoes within the rectangle:
volcidstoreplace = list(v.get_volcanoes_in_polygon(newpoly.geom[0]).volc_id.values)
print(volcidstoreplace)
cmds = []
print('printing only for a check:')
for secvolc in volcidstoreplace:
    if not secvolc == primvolc:
        oldvclip = v.get_volclip_vids(secvolc)[0]
        oldpath = os.path.join(os.environ['LiCSAR_volc'], str(oldvclip))
        if os.path.exists(oldpath):
            #cmd = 'rm -rf '+oldpath
            cmd = "mv $LiCSAR_volc/{0} $LiCSAR_volc/../volc_subsets.backup/{1}".format(str(oldvclip), str(secvolc))
            if oldvclip > 2700:
                print('Not doing :'+str(cmd))
            else:
                print(cmd)
                cmds.append(cmd)


# if all ok, then:
if cmds:
    for c in cmds:
        print(c)
        os.system(c)


for secvolc in volcidstoreplace:
    if not secvolc == primvolc:
        try:
            oldvclip = v.get_volclip_vids(secvolc)[0]
            if oldvclip > 2700:
                print('Not removing old vclip for volcano '+str(secvolc))
            else:
                v.delete_volclip(oldvclip)
        except:
            print('no old clip for volcano id '+str(secvolc))
        v.add_volcano_to_volclip(secvolc, newclip)


### 4. initialise new clip:

checkfrs = v.get_volcano_frames(vid)
print('import volcdb as v')
for chfr in checkfrs:
    print('v.initialise_subset_volclip(vid = '+str(newclip)+', frame = "'+chfr+'")')
    # v.initialise_subset_volclip(vid=newclip, frame=chfr)




v.initialise_subset_volclip(vid = 2830, frame = "060A_00001_030604")
v.initialise_subset_volclip(vid = 2830, frame = "169D_00001_020800")

v.initialise_subset_volclip(vid = 2833, frame = "013D_04281_131313")
v.initialise_subset_volclip(vid = 2833, frame = "064A_04218_131313")
v.initialise_subset_volclip(vid = 2833, frame = "115D_04405_121313")
v.initialise_subset_volclip(vid = 2833, frame = "137A_04234_141413")
