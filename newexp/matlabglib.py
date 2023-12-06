import scipy.io
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from tqdm import tqdm
from scipy.ndimage import label 
from scipy.ndimage import binary_dilation as bd
import copy
import sys

def generateBedmapGLIBs(depth,ice,resolution=5):
    #We don't care about the structure of the bathymetry off the shelf break for now
    # so we'll call the off shelf water the part that's greater than 2000 meters deep
    # and just set it to 2000 meters
    lowest= (-np.nanmin(depth))-100
    depth[depth<-lowest] = -lowest

    ## Grab a single point in the open ocean. We just made all the open ocean one depth so it is connected 
    openoceancoord = np.where(depth==-lowest)
    openoceancoord = (openoceancoord[0][0],openoceancoord[1][0])

    ## Decide the vertical resolution at which we'd like to search from -2000 to 0
    slices = - np.asarray(range(0,int(lowest),resolution)[::-1])
    ## A matrix output we will fill with nans as we go.
    GLIB = np.full_like(depth,np.nan)
    for iso_z in slices:
        ### A mask that shows all points below our search depth that arent land or grounded ice.
        below = np.logical_and((depth<=iso_z),ice!=0)
        ### Use the image library label function to label connected regions.
        regions, _ = label(below)
        ### If any region is connected to open ocean at this depth and it hasn't been assigned a GLIB before (all ocean points are connected at z=0)
        ### then set it's glib to iso_z -20. The minus 20 is there because we are looking for the depth at which its not connected
        GLIB[np.logical_and(regions==regions[openoceancoord],np.isnan(GLIB))] = iso_z-resolution
    return GLIB


def GLIBfromFile(fname,full=False):
    variables = scipy.io.loadmat(fname,variable_names=('h','icedraft'))
    icedraft = np.asarray(variables["icedraft"])
    h = np.asarray(variables["h"])
    icedraft[icedraft==0] = np.nan
    icedraft[h==0] = 0
    GLIB = generateBedmapGLIBs(h,icedraft)
    GLIB[np.asarray(GLIB-h)<10] = np.nan
    #plt.imshow(GLIB)
    #plt.show()
    if full:
        return GLIB
    else:
        return np.nanmean(GLIB[icedraft==np.nanmin(icedraft)])

print(GLIBfromFile(sys.argv[1]))
