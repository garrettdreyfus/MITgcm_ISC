import scipy.io
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from tqdm import tqdm
from scipy.ndimage import label 
from scipy.ndimage import binary_dilation as bd


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
    for iso_z in tqdm(slices):
        ### A mask that shows all points below our search depth that arent land or grounded ice.
        below = np.logical_and((depth<=iso_z),ice!=0)
        ### Use the image library label function to label connected regions.
        regions, _ = label(below)
        ### If any region is connected to open ocean at this depth and it hasn't been assigned a GLIB before (all ocean points are connected at z=0)
        ### then set it's glib to iso_z -20. The minus 20 is there because we are looking for the depth at which its not connected
        GLIB[np.logical_and(regions==regions[openoceancoord],np.isnan(GLIB))] = iso_z-resolution
    return GLIB


def calcMSD(depth,ice,slice_depth,resolution=5):
    #We don't care about the structure of the bathymetry off the shelf break for now
    # so we'll call the off shelf water the part that's greater than 2000 meters deep
    # and just set it to 2000 meters
    lowest= (-np.nanmin(depth))-100
    depth[depth<-lowest] = -lowest

    ## Grab a single point in the open ocean. We just made all the open ocean one depth so it is connected 
    openoceancoord = np.where(depth==-lowest)
    openoceancoord = (openoceancoord[0][0],openoceancoord[1][0])

    ## Decide the vertical resolution at which we'd like to search from -2000 to 0
    ## A matrix output we will fill with nans as we go.
    MSD = np.full_like(depth,np.nan)
    below = ~np.logical_and((depth<=slice_depth),ice!=0)
    regions, _ = label(~below)
    for iters in range(0,100,resolution):
        ### A mask that shows all points below our search depth that arent land or grounded ice.
        ### Use the image library label function to label connected regions.
        regions, _ = label(~below)
        ### If any region is connected to open ocean at this depth and it hasn't been assigned a GLIB before (all ocean points are connected at z=0)
        ### then set it's glib to iso_z -20. The minus 20 is there because we are looking for the depth at which its not connected
        condition_part1 = np.logical_and(regions!=regions[openoceancoord],np.isnan(MSD))
        condition = np.logical_and(condition_part1,below!=1)
        MSD[condition] = iters
        fig, (ax1,ax2,ax3) = plt.subplots(1,3)
        ax1.imshow(below)
        ax2.imshow(regions)
        ax3.imshow(MSD)
        plt.show()
        below = bd(below,iterations=resolution,mask=regions)
    #plt.imshow(MSD)
    #plt.show()
    return MSD

#variables = scipy.io.loadmat('../experiments/GLIB-explore/under/input/metaparameters.mat',variable_names=('h','icedraft'))
variables = scipy.io.loadmat('../experiments/smallerslope-GLIB-explore-18/at-125/input/metaparameters.mat',variable_names=('h','icedraft'))
icedraft = np.asarray(variables["icedraft"])
h = np.asarray(variables["h"])
plt.imshow(h)
plt.show()
icedraft[icedraft==0] = np.nan
icedraft[h==0] = 0

plt.imshow(icedraft==np.nanmin(icedraft))
plt.show()


GLIB = generateBedmapGLIBs(h,icedraft)
#MSD = calcMSD(h,icedraft,-500)
#MSD = calcMSD(h,icedraft,-200)

randomcmap = matplotlib.colors.ListedColormap(np.random.rand ( 256,3))
GLIB[np.asarray(GLIB-h)<10] = np.nan
print(np.nanmean(GLIB[icedraft==np.nanmin(icedraft)]))
plt.imshow(GLIB,cmap="jet")#,cmap=randomcmap)
plt.colorbar()

plt.show()

