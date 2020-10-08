#!/usr/bin/env python

import sys
import numpy as np
from PIL import Image
#import time
#from matplotlib import pyplot as plt
import cv2
import os
#import glob
from skimage import morphology
from osgeo import osr, gdal
import datetime
import subprocess as subp


#try:
#    wrap=sys.argv[1]
#    unwrap=sys.argv[2]
#except:
#    print('please provide parameters: wrap_ifg_file unwrap_ifg_file')
#    exit()

def qc1(a,image,bound):
    thresh = cv2.threshold(a, 10, 1, cv2.THRESH_BINARY)[1]  #  1 was 255
    kernel=cv2.getStructuringElement(cv2.MORPH_RECT,(15,15))
    thresh = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, kernel)
    thresh=thresh.astype(bool)
    thresh = morphology.remove_small_holes(thresh, 200)
    thresh = morphology.remove_small_objects(thresh, 200)
    thresh=np.uint8(1*thresh) # convert bol into int # was 255 instead of 1
    ###### line detection
    low_threshold = 0.2 # was 50
    high_threshold = 0.58 # was 150
    kernel=5
    edges = cv2.Canny(thresh, low_threshold, high_threshold,kernel)
    edges = np.multiply(bound, edges)
    # plt.figure()
    # plt.imshow(edges, cmap='gray')
    ####### Hough line detection
    rho = 1  # distance resolution in pixels of the Hough grid
    theta = np.pi / 180  # angular resolution in radians of the Hough grid
    threshold = 10 # minimum number of votes (intersections in Hough grid cell)
    min_line_length = 80  # minimum number of pixels making up a line
    max_line_gap = 5  # maximum gap in pixels between connectable line segments
    line_image = np.copy(image) * 0  # creating a blank to draw lines on
    # Run Hough on edge detected image
    # Output "lines" is an array containing endpoints of detected line segments
    lines = cv2.HoughLinesP(edges, rho, theta, threshold, np.array([]),
                        min_line_length, max_line_gap)
    if lines is None:
        return(0)
    return (len(lines))


def check_dimensions(tiffile1, tiffile2):
    #flag 0 means ok
    try:
        rds1 = gdal.Open(tiffile1)
        rds2 = gdal.Open(tiffile2)
    except:
        print('error reading file')
        return 1
    if (rds1.RasterXSize == rds2.RasterXSize and rds1.RasterYSize == rds2.RasterYSize):
        flag = 0
    else:
        flag = 1
    return flag


def check_lines(png):
    img = Image.open(png)
    data = np.array( img, dtype='uint8')
    image = cv2.imread(png)
    aa= cv2.cvtColor(image,cv2.COLOR_BGR2GRAY)
    try:
        a=data[:,:,3]
    except:
        a = aa
    a_inv=255-a
    #imagew = cv2.imread(wrap)
    thresh = cv2.threshold(aa, 0, 255, cv2.THRESH_BINARY)[1]
    kernel=cv2.getStructuringElement(cv2.MORPH_RECT,(5,5))
    mask = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, kernel)
    mask_dil = cv2.morphologyEx(mask, cv2.MORPH_DILATE, kernel)
    mask_erd = cv2.morphologyEx(mask, cv2.MORPH_ERODE, kernel)
    bound = mask_dil - mask_erd
    bound = cv2.morphologyEx(bound, cv2.MORPH_DILATE, kernel)
    bound=np.uint8(bound/255)
    bound = 1 - bound
    
    im1=qc1(a,image,bound) # better for boundary extraction
    im2=qc1(a_inv,image,bound)
    n_lines=im1+im2
    if (n_lines>3):
        flag = 1
    else:
        flag = 0
    return flag


def check_lines_orig(wrap, unwrap):
    #global wrap
    #global unwrap
    img = Image.open(unwrap)
    data = np.array( img, dtype='uint8')
    image = cv2.imread(unwrap)
    a=data[:,:,3]
    a_inv=255-a
    imagew = cv2.imread(wrap)
    aa= cv2.cvtColor(imagew,cv2.COLOR_BGR2GRAY)
    thresh = cv2.threshold(aa, 0, 255, cv2.THRESH_BINARY)[1]
    kernel=cv2.getStructuringElement(cv2.MORPH_RECT,(5,5))
    mask = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, kernel)
    mask_dil = cv2.morphologyEx(mask, cv2.MORPH_DILATE, kernel)
    mask_erd = cv2.morphologyEx(mask, cv2.MORPH_ERODE, kernel)
    bound = mask_dil - mask_erd
    bound = cv2.morphologyEx(bound, cv2.MORPH_DILATE, kernel)
    bound=np.uint8(bound/255)
    bound = 1 - bound
    
    im1=qc1(a,image,bound) # better for boundary extraction
    im2=qc1(a_inv,image,bound)
    n_lines=im1+im2
    if (n_lines>3):
        flag = 1
    else:
        flag = 0
    return flag


def check_timescan(check_coh_file, thres1 = 0.01):
    badifgs_timescan = []
    #now check using the full series:
    if os.path.exists(check_coh_file):
        data = np.genfromtxt(check_coh_file, unpack=True).T
        filedate = np.genfromtxt(check_coh_file, usecols=0, dtype=str)
        Bt=data[:,1]
        coh_mean=data[:,4]
        frac_unw=data[:,3]
        frac_data=data[:,2]
        max_coh=np.max(coh_mean)
        max_frac=np.max(frac_unw)
        max_data=np.max(frac_data)
        Q0=np.round(100*frac_data/max_data)
        Q1=coh_mean/max_coh
        Q2=frac_unw/max_frac
        Btmin=np.min(Bt)
        Btmax=5*Btmin
        #this is to identify bad processed ifg
        #thres1 = coh_threshold
        bad_ind=np.logical_or(Q1<thres1,Q2<thres1)
        bad_processed = filedate[bad_ind]
        data1_save = np.array([bad_processed])
        data1_save = data1_save.T
        for ifg in data1_save:
            if len(ifg) == 1:
                ifg = ifg[0]
            badifgs_timescan.append(ifg)
    else:
        print('error - the coh stat file {} does not exist'.format(check_coh_file))
        return []
    return badifgs_timescan

def duration(tmaster, tslave):
    t1 = datetime.datetime.strptime(tmaster, "%Y%m%d")
    t2 = datetime.datetime.strptime(tslave, "%Y%m%d")
    return (t2 - t1).days


def func1(lim1, lim2):
    ind=np.logical_and(Q0>lim1,Q0<=lim2)
    q1 = Q1[ind]
    q2 = Q2[ind]
    bt = Bt[ind]
    return q1, q2, bt
 

def get_stats(path, ifg):
    #cctif, unwtif, ):
    tmaster = ifg[0:8]
    tslave =  ifg[9:17]
    Bt = duration(tmaster,tslave)
    cctif = os.path.join(path,ifg+'.geo.cc.tif')
    unwtif = os.path.join(path,ifg+'.geo.unw.tif')
    try:
        gtif = gdal.Open(cctif)
    except IOError:
        #print('The file {} is corrupted'.format(cctif))
        return False
    try:
        test = gdal.Open(unwtif)
    except IOError:
        #print('The file {} is corrupted'.format(unwtif))
        return False
    if test==None:
        #corrupted_file
        return False
    data = gtif.ReadAsArray()
    Ntotal = np.size(data)
    if Ntotal<2:
        Frac=0
    else:
        num = np.where((abs(data)>0.0))
        Nproc = np.size(num)/2
        Frac = Nproc/Ntotal
    Mean = np.mean(data)
    Std = np.std(data)
    Min = np.amin(data)
    Max = np.amax(data)
    
    gtif2 = gdal.Open(unwtif)
    data2 = gtif2.ReadAsArray()
    Ntotal = np.size(data2)
    if Ntotal<2:
        Frac2=0
    else:
        num2 = np.where((abs(data2)>0.0))
        Nunw = np.size(num2)/2
        Frac2 = Nunw/Ntotal
    #stats = ""
    stats = "%s %d %.3f %.3f %.3f %.3f %.3f %.3f\n" % (ifg,Bt,Frac,Frac2,Mean,Std,Min,Max)
    return stats


def basic_check(ifgdir):
    #returns True for 'ifg files are ok'
    output = True
    #types = ['cc','di']'
    needed_files_ext = [ 'cc.png', 'cc.tif', 'diff_pha.tif', 'diff.png', 'unw.png', 'unw.tif']
    ifg = os.path.basename(ifgdir)
    for ext in needed_files_ext:
        if not os.path.exists(os.path.join(ifgdir,ifg+'.geo.'+ext)):
            return False
        else:
            a = subp.run(['gdalinfo', os.path.join(ifgdir,ifg+'.geo.'+ext)], capture_output=True)
            if 'ERROR' in a.stderr.decode():
                return False
            elif os.path.getsize(os.path.join(ifgdir,ifg+'.geo.'+ext)) == 0:
                return False
    return output


def main():
    try:
        ifgdir = sys.argv[1]
    except:
        print('please provide parameter: interferograms directory, e.g. ')
        print('quality_check.py /gws/nopw/j04/nceo_geohazards_vol1/public/LiCSAR_products/144/144D_05102_121313/interferograms')
        exit()
    if not os.path.exists(ifgdir):
        print('the provided directory {} does not exist'.format(ifgdir))
        exit()
    print('this is using both line detection and timescan approaches')
    print('bad interferograms: ')
    badifgs = []
    check_coh_file = os.path.join(ifgdir, 'check_coherence.txt')
    if os.path.exists(check_coh_file):
        os.remove(check_coh_file)
    coh_threshold = 0.01
    
    
    for ifg in os.listdir(ifgdir):
        #check the lines in wrapped imgs
        #check only wrapped imgs:
        wrap = os.path.join(ifgdir,ifg,ifg+'.geo.diff.png')
        unwrap = os.path.join(ifgdir,ifg,ifg+'.geo.unw.png')
        cctif = os.path.join(ifgdir,ifg,ifg+'.geo.cc.tif')
        if not basic_check(os.path.join(ifgdir,ifg)):
            #if (not os.path.exists(wrap)) or (not os.path.exists(unwrap)) or (not os.path.exists(cctif)):
            flag = 1
        else:
            flag = check_lines(wrap) #, unwrap)
        
        if flag == 0:
            #check 
            stats = get_stats(os.path.join(ifgdir,ifg), ifg) 
            if not stats:
                flag = 1
            else:
                fid = open(check_coh_file, 'a')
                fid.write(stats)
                fid.close()
        if flag == 1:
            badifgs.append(ifg)
    #just print the bad ifgs now
    for ifg in badifgs:
        print(ifg)
    badifgs = check_timescan(check_coh_file, coh_threshold)
    for ifg in badifgs:
        print(ifg)
    if os.path.exists(check_coh_file):
        os.remove(check_coh_file)

if __name__ == '__main__':
    main()
