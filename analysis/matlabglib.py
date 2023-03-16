import scipy.io
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from tqdm import tqdm
from scipy.ndimage import label 
from scipy.ndimage import binary_dilation as bd
import copy

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

def calcMSD(depth,ice,slice_depth,points_of_interest,resolution=5,debug=False):

    # First let set the open ocean all to one depth
    lowest= (-np.nanmin(depth))-100
    depth[depth<-lowest] = -lowest
    depth[np.isnan(depth)] = -lowest

    ## Grab a single point in the open ocean. We just made all the open ocean one depth so it is connected 
    openoceancoord = ([145],[100])

    ## We're going to need to move the points of interest
    ## this is the drawstring bag approach of course
    moving_poi=copy.deepcopy(points_of_interest)

    ## Here's where we'll store actual results
    MSD = points_of_interest.shape[0]*[np.nan]

    ## And let's make sure to count how many points are even eligible for finding the MSD
    ## They are eligible if the points are actually above bed rock at this depth.
    count = 0
    valid_poi = []
    for l in range(len(MSD)):
        if depth[points_of_interest[l][0],points_of_interest[l][1]]<slice_depth:
            count+=1
            valid_poi.append(True)
        else:
            valid_poi.append(False)

    valid_poi = np.asarray(valid_poi)
    ## A binary mask of all points below the slice
    below = ~np.logical_and((depth<=slice_depth),ice!=0)
    border_mask = np.full_like(depth,1)
    border_mask[:,-10:] = 0
    dilation = 0


    while np.sum(valid_poi)>0:
        # print(np.sum(valid_poi))

        ## we now calculate the connectedness
        regions, _ = label(~below)

        ## We need to move the points of interest out of
        ## bedrock if they have been dilated within

        for idx in range(len(moving_poi)):
            if regions[moving_poi[idx][0],moving_poi[idx][1]]==0 and valid_poi[idx]:
                i_idx = moving_poi[idx][0]
                j_idx = moving_poi[idx][1]
                r = resolution+1
                left = max(i_idx-r,0)
                right = min(i_idx+r+1,depth.shape[0])
                down = max(j_idx-r,0)
                up = min(j_idx+r+1,depth.shape[1])
                neighborhood = regions[left:right,down:up]
                if np.max(neighborhood)>0:
                    offset = np.where(neighborhood!=0)
                    offset_i = offset[0][0]-(i_idx-left)
                    offset_j = offset[1][0]-(j_idx-down)
                    before = regions[moving_poi[idx][0],moving_poi[idx][1]]
                    moving_poi[idx][0] = moving_poi[idx][0]+offset_i
                    moving_poi[idx][1] = moving_poi[idx][1]+offset_j
                    after = regions[moving_poi[idx][0],moving_poi[idx][1]]
                else:
                    print(neighborhood)
                    print(regions[left-5:right+5,down-5:up+5])
                    MSD[idx]=dilation
                    valid_poi[idx]=False
                    print("zonk")
                        
        ## now we check for connectedness
        for idx in range(len(moving_poi)):
            if regions[moving_poi[idx][0],moving_poi[idx][1]] != regions[openoceancoord[1][0],openoceancoord[0][0]] and \
               valid_poi[idx]:
                MSD[idx]=dilation
                valid_poi[idx]=False

        if debug:
            plt.imshow(regions)
            plt.scatter(moving_poi.T[1],moving_poi.T[0],c="red")
            print(openoceancoord)
            plt.scatter(openoceancoord[0],openoceancoord[1],c="red")
            plt.colorbar()
            plt.show()
        
        below = bd(below,iterations=resolution,border_value=1,mask=border_mask)
        dilation+=resolution

    if debug:
        print(MSD)
    # plt.imshow(np.logical_and((depth<=slice_depth),ice!=0))
    # plt.scatter(points_of_interest.T[1],points_of_interest.T[0],c=MSD,cmap="jet")
    # plt.colorbar()
    # plt.show()
    return MSD

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

def aisfdepth(fname,full=False):
    variables = scipy.io.loadmat(fname,variable_names=('h','icedraft'))
    icedraft = np.asarray(variables["icedraft"])
    h = np.asarray(variables["h"])
    icedraftbd = bd(icedraft!=0)
    #fig, (ax1,ax2) = plt.subplots(1,2)
    #ax1.imshow(icedraft)
    #ax2.imshow(np.logical_and(icedraftbd,np.logical_and(icedraft==0,h!=0)))
    #ax2.imshow(h)
    #plt.show()
    return np.nanmean(h[np.logical_and(icedraftbd,np.logical_and(icedraft==0,h!=0))])

def MSDfromFile(fname):
    print(fname)
    variables = scipy.io.loadmat(fname,variable_names=('h','icedraft','dx','dy'))
    icedraft = np.asarray(variables["icedraft"])
    dx,dy = np.asarray(variables["dx"])[0],np.asarray(variables["dy"])[0]
    dist = np.mean([dx,dy])
    h = np.asarray(variables["h"])
    GLIB=GLIBfromFile(fname,True)
    icedraft[icedraft==0] = np.nan
    icedraft[h==0] = 0
    poi = np.asarray(np.where(icedraft==np.nanmin(icedraft))).T
    x,y = np.where(icedraft==np.nanmin(icedraft))
    start = int(np.nanmin(GLIB[x,y]))
    depths = range(start,0,20)
    MSDs = []
    for k in depths:
        MSD = calcMSD(h,icedraft,k,poi,resolution=1,debug=False)#(k==-200))
        MSDs.append(np.nanmedian(MSD)*dist/1000)
    return MSDs,depths

def RdfromFile(fname):
    variables = scipy.io.loadmat(fname,variable_names=('Rd','yy'))
    Rd = np.asarray(variables["Rd"])
    yy = np.asarray(variables["yy"])


#m,d = MSDfromFile('../experiments/halfw-GLIB-explore-18/at-125/input/metaparameters.mat')
#RdfromFile('../experiments/GLIB-explore-18/at-125/input/metaparameters.mat')
#plt.plot(m,d,c="orange")
#m,d = MSDfromFile('../experiments/doublew-GLIB-explore-18/at-125/input/metaparameters.mat')
#plt.plot(m,d,c="gray")
#m,d = MSDfromFile('../experiments/GLIB-explore-16/at-125/input/metaparameters.mat')
#plt.plot(m,d,c="red")
#m,d = MSDfromFile('../experiments/GLIB-explore-18/at-125/input/metaparameters.mat')
#plt.plot(m,d,c="green")
#plt.xlabel("MSD")
#plt.ylabel("Depth")
#plt.show()


#print("700",GLIBfromFile('../experiments/shelfzexp-GLIB-explore-101/at0d700/input/metaparameters.mat'))
#print("500",GLIBfromFile('../experiments/shelfzexp-GLIB-explore-101/at0d500/input/metaparameters.mat'))

#variables = scipy.io.loadmat('../experiments/GLIB-explore/under/input/metaparameters.mat',variable_names=('h','icedraft'))
# variables = scipy.io.loadmat('../experiments/smallerslope-GLIB-explore-18/at-125/input/metaparameters.mat',variable_names=('h','icedraft'))
# icedraft = np.asarray(variables["icedraft"])
# h = np.asarray(variables["h"])
# plt.imshow(h)
# plt.show()
# icedraft[icedraft==0] = np.nan
# icedraft[h==0] = 0

# plt.imshow(icedraft==np.nanmin(icedraft))
# plt.show()


# GLIB = generateBedmapGLIBs(h,icedraft)
# #MSD = calcMSD(h,icedraft,-500)
# #MSD = calcMSD(h,icedraft,-200)

# randomcmap = matplotlib.colors.ListedColormap(np.random.rand ( 256,3))
# GLIB[np.asarray(GLIB-h)<10] = np.nan
# print(np.nanmean(GLIB[icedraft==np.nanmin(icedraft)]))
# plt.imshow(GLIB,cmap="jet")#,cmap=randomcmap)
# plt.colorbar()

# plt.show()

