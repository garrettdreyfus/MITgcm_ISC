import numpy as np
import glob
from xmitgcm import open_mdsdataset
import numpy.matlib 
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegFileWriter
from tqdm import tqdm



def timeSeriesDashboard(fname,label,fig,axises,times=np.array([])):
    ((ax1,ax2,ax5),(ax3,ax4,ax6)) = axises 
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    if times.any():
        ds = open_mdsdataset(fname,ignore_unknown_vars=True,iters=times,extra_variables = extra_variables)
    else:
        times=getIterNums(fname)
        print(times)
        ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    totalvolume  = (ds.XC*ds.YC*ds.Z*ds.hFacC).sum().values
    ## theta plot

    ts = ds.time.values*108.0/60.0/60.0/24.0/365.0
    tsnew = np.full_like(ts,0,dtype=float)
    tsnew[:] = ts
    ts = tsnew
    ts = ts/1000000000
    times=ts
    totaltheta = (ds.THETA*ds.XC*ds.YC*ds.Z*ds.hFacC).sum(axis=[1,2,3]).values
    ax1.plot(times,totaltheta/totalvolume,label=label)
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Potential Temperature")
    ax1.legend()
    ## salt plot
    totalsalt = (ds.SALT*ds.XC*ds.YC*ds.Z*ds.hFacC).sum(axis=[1,2,3]).values
    ax2.plot(times,totalsalt/totalvolume,label=label)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Salinity")

    ## kinetic energy plot
    totalmomke = (ds.momKE*ds.XC*ds.YC*ds.Z*ds.hFacC).sum(axis=[1,2,3]).values
    ax3.plot(times,totalmomke/totalvolume,label=label)
    ax3.set_xlabel("Time")
    ax3.set_ylabel("Kinetic Energy")

    
    ## kinetic energy plot
    #totalSHIfwFlx = (ds.SHIfwFlx*ds.XC*ds.YC*ds.Z*ds.hFacC).sum(axis=[1,2,3]).values
    mask = ~np.isnan(ds.SHIfwFlx.values)
    shflx = ds.SHIfwFlx.values
    print(shflx)
    totalSHIfwFlx = np.sum(shflx*mask,axis=(1,2))*(60*60*24*365)*(1/920.0)*(1/np.sum(mask,axis=(1,2)))
    print(ds.SHIfwFlx)
    ax4.plot(times,totalSHIfwFlx,label=label)
    ax4.set_xlabel("Time")
    ax4.set_ylabel("Fresh Water Flux")
    ax6.remove()

    bottomtemps = []
    surfacetemps = []
    for k in range(ds["THETA"].shape[0]):
        z = np.concatenate((ds.hFacC.values[:-1,:,:]-ds.hFacC.values[1:,:,:],ds.hFacC[-1:,:,:]),axis=0)
        z[z<0]=0
        z[-1,:,:]=0
        zmask = z
        X = np.full_like(ds["THETA"][k].values,np.nan,dtype=float)
        X[ds.hFacC.values != 0]= ds["THETA"][k].values[ds.hFacC.values != 0]
        znew = np.multiply(zmask,X)
        nancount = np.nansum(np.isnan(znew),axis=0)
        znew = np.nansum(znew,axis=0)
        znew[nancount==X.shape[0]] = np.nan
        bottomtemps.append(np.nanmean(znew))
        surfaceslice = ds["THETA"][k][0,:,:].values
        surfaceslice[ds.hFacC[0,:,:].values==0]=np.nan
        surfaceslice = np.nanmean(surfaceslice)
        surfacetemps.append(np.nanmean(surfaceslice))
    ax5.plot(times,bottomtemps,label="Bottom")
    ax5.plot(times,surfacetemps,label="Surface")
    ax5.set_xlabel("Time")
    ax5.set_ylabel("Potential Temperature")
    ax5.legend()

def depthFromdZ(ds):
    U = ds.UVEL.values[0,:,:,:]
    fZ = list(ds.Z)
    DZ = np.asarray(fZ)
    Z = np.repeat(DZ,U.shape[1]*U.shape[2])
    Z = Z.reshape(U.shape[0],U.shape[1],U.shape[2])
    z = np.concatenate((ds.hFacC.values[:-1,:,:]-ds.hFacC.values[1:,:,:],ds.hFacC[-1:,:,:]),axis=0)
    z[z<0]=0
    z[-1,:,:]=0
    zmask = z
    return np.sum(zmask*Z,axis=0)


def barotropic_streamfunction(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    shortname, fpath = outPath(fname) 
    moviewriter = FFMpegFileWriter(fps=1)
    fig,ax1 = plt.subplots()
    depth = depthFromdZ(ds)
    UFULL = ds.UVEL.values
    VFULL = ds.VVEL.values
    with moviewriter.saving(fig, fpath+"-bt.mp4" , dpi=250):
        for k in tqdm(range(0,ds.UVEL.values.shape[0])):
            U = UFULL[k,:,:,:]
            V = VFULL[k,:,:,:]
            fDZ = list(np.diff(ds.Z))
            fDZ.append(fDZ[-1])
            fDZ = np.asarray(fDZ)
            DZ = np.repeat(fDZ,U.shape[1]*U.shape[2])
            DZ = DZ.reshape(U.shape[0],U.shape[1],U.shape[2])
            UU = np.sum(U*DZ*ds.hFacW.values,axis=0)
            xs = ds.XC.values
            ys = list(np.diff(ds.YC.values))
            ys.append(ys[-1])
            ys = np.repeat(ys,U.shape[2])
            ys = ys.reshape(U.shape[1],U.shape[2],order="F")
            bt = np.cumsum(UU*ys,axis=0)
            bt[np.sum(ds.hFacC,axis=0)==0] = np.nan
            frame = ax1.pcolormesh(ds.XC.values,ds.YC.values,bt)
            ax1.contour(ds.XC.values,ds.YC.values,depth,levels=20,colors="black")
            #ax1.quiver(ds.XC.values,ds.YC.values,np.sum(U,axis=0),np.sum(V,axis=0))
            plt.title(str(k))
            cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            cb.remove()
            frame.remove()

def meltmap(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    shortname, fpath = outPath(fname) 
    moviewriter = FFMpegFileWriter(fps=1)
    fig,ax1 = plt.subplots()
    print(ds)
    mask = np.logical_and(ds.hFacC.values[0]==0,np.sum(ds.hFacC.values,axis=0)!=0)
    with moviewriter.saving(fig, fpath+"-meltmap.mp4" , dpi=250):
        for k in tqdm(range(ds.UVEL.values.shape[0])):
            melt= ds.SHIfwFlx.values[k]
            melt[~mask]=np.nan
            frame = ax1.pcolormesh(melt)
            cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            cb.remove()
            frame.remove()

def bottomMask(ds):
    bmask = np.full_like(ds.THETA[0,0,:,:],0,dtype=int)
    print(ds)
    for k in range(ds.Z.shape[0])[::-1]:
        nanmask = (bmask == 0)
        full_mask = np.full_like(bmask,k,dtype=int)
        full_mask[ds.maskC[k,:,:]==0] = 0
        bmask[nanmask] = full_mask[nanmask]
    return bmask
    
def outPath(resultspath):
    nameparts = resultspath.split("/")
    shortname = nameparts[-3] + "|" + nameparts[-2]
    fpath = "/".join((nameparts[:-4]+["pics"] + [shortname]))
    return shortname, fpath

def crossSectionAnim(fname,description,times=np.array([]),quant="THETA"):
    fig,ax1 = plt.subplots()
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    print(times)
    ds = open_mdsdataset(fname,prefix=quant,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    ds[quant].values=ds[quant].values*ds.hFacC.values
    #print(ds.hFacC)
    zonal_average = ds.where(ds.hFacC ==1).mean(dim="XC",skipna=True)
    #zonal_average = ds.mean(dim="XC",skipna=True)
    #zonal_average = ds.isel(XC=50)
    moviewriter = FFMpegFileWriter(fps=1)
    tmin, tmax = np.nanmin(zonal_average[quant]), np.nanmax(zonal_average[quant])
    shortname, fpath = outPath(fname) 
    print(shortname,fpath)
    fig.suptitle(shortname)
    ys = zonal_average.YC.values
    zs = zonal_average.Z.values
    zvals = zonal_average[quant].values
    with moviewriter.saving(fig, fpath+".mp4" , dpi=250):
        print("writing movie")
        for k in tqdm(range(zonal_average[quant].shape[0])):
            frame = ax1.pcolormesh(ys,zs,zvals[k,:,:],cmap="jet",vmin=tmin,vmax=tmax)
            cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            cb.remove()
            frame.remove()

def bottomAnim(fname,description,times=np.array([]),quant="SALT"):
    fig,ax1 = plt.subplots()
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    z = np.concatenate((ds.hFacC.values[:-1,:,:]-ds.hFacC.values[1:,:,:],ds.hFacC[-1:,:,:]),axis=0)
    z[z<0]=0
    z[-1,:,:]=0
    zmask = z
    #zonal_average = ds.isel(XC=32)
    moviewriter = FFMpegFileWriter(fps=1)
    tmin, tmax = np.nanmin(ds[quant]), np.nanmax(ds[quant])
    shortname, fpath = outPath(fname) 
    depth = depthFromdZ(ds)
    quantvals = ds[quant].values
    with moviewriter.saving(fig, fpath+"-bot.mp4", dpi=250):
        print("writing movie")
        for k in tqdm([0]+list(range(quantvals.shape[0]))+[-1]):
            d = np.mean(quantvals,axis=0)
            X = np.full_like(d,np.nan,dtype=float)
            X[ds.hFacC.values != 0]= d[ds.hFacC.values != 0]
            znew = np.multiply(zmask,X)
            nancount = np.nansum(np.isnan(znew),axis=0)
            znew = np.nansum(znew,axis=0)
            znew[nancount==X.shape[0]] = np.nan
            frame = ax1.pcolormesh(ds.XC.values,ds.YC.values,znew,cmap="jet",vmin=34.4,vmax=34.65)
            ax1.contour(ds.XC.values,ds.YC.values,depth,colors="black",levels=20)
            cb = plt.colorbar(frame)
            plt.show()
            moviewriter.grab_frame()
            cb.remove()
            frame.remove()


def surfaceAnim(fname,description,times=np.array([]),quant="SALT"):
    fig,ax1 = plt.subplots()
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    z = np.concatenate((ds.hFacC.values[:-1,:,:]-ds.hFacC.values[1:,:,:],ds.hFacC[-1:,:,:]),axis=0)
    z[z>0]=0
    z[-1,:,:]=0
    zmask = z
    #zonal_average = ds.isel(XC=32)
    moviewriter = FFMpegFileWriter(fps=1)
    tmin, tmax = np.nanmin(ds[quant]), np.nanmax(ds[quant])
    shortname, fpath = outPath(fname) 
    depth = depthFromdZ(ds)
    quantvals = ds[quant].values
    with moviewriter.saving(fig, fpath+"-surf.mp4", dpi=250):
        print("writing movie")
        for k in tqdm([0]+list(range(0,quantvals.shape[0]))+[-1]):
            X = np.full_like(quantvals[k],np.nan,dtype=float)
            X[ds.hFacC.values != 0]= quantvals[k][ds.hFacC.values != 0]
            znew = np.multiply(zmask,X)
            nancount = np.nansum(np.isnan(znew),axis=0)
            znew = np.nansum(znew,axis=0)
            znew[nancount==X.shape[0]] = np.nan
            frame = ax1.pcolormesh(ds.XC.values,ds.YC.values,znew,cmap="jet",vmin=-2,vmax=1)
            ax1.contour(ds.XC.values,ds.YC.values,depth,colors="black",levels=20)
            cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            cb.remove()
            frame.remove()



def getIterNums(fpath):
    iters = []
    saltiters = []
    for fname in glob.glob(fpath+"/*.data"):
        if 'THETA' in fname and "inst" not in fname:
            n = int(fname.split(".")[1])
            iters.append(n)
        if 'SALT' in fname and "inst" not in fname:
            n = int(fname.split(".")[1])
            saltiters.append(n)
    return np.unique(np.asarray(np.intersect1d(iters,saltiters)))[:-1]

#timeSeriesDashboard("/home/garrett/Projects/MITgcm_ISC/experiments/squish/test/results","Restricted y domain length with default settings",times = np.asarray(range(1,9))*420480)
#fig,axises = plt.subplots(2,3)
#timeSeriesDashboard("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-2O-rand/results","2O-r",fig,axises)
#timeSeriesDashboard("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG/results","PIG",fig,axises)
#timeSeriesDashboard("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-2O/results","2O",fig,axises)
#timeSeriesDashboard("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-steep/results","steep",fig,axises)
#plt.show()

#meltmap("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG/results","Restricted y domain length with default settings")
#barotropic_streamfunction("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-steep/results","Restricted y domain length with default settings")
#barotropic_streamfunction("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-2O/results","Restricted y domain length with default settings")
#barotropic_streamfunction("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-2O-rand/results","Restricted y domain length with default settings")
#barotropic_streamfunction("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-2O-rand/results","Restricted y domain length with default settings")
#barotropic_streamfunction("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG/results","Restricted y domain length with default settings")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/squish-polyna/under/results","Restricted y domain length with default settings")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/fully/under/results","Restricted y domain length with default settings")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/orlanski-test/at/results","Restricted y domain length with default settings")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/squish-polyna/above/results","Restricted y domain length with default settings",quant="SALT")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/squish-polyna/at/results","Restricted y domain length with default settings",quant="SALT")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/randtopo/crashtest1/results","Restricted y domain length with default settings",quant="THETA")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/randtopo/crashtest-orlanskiw/results","Restricted y domain length with default settings",quant="THETA")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/orlanski-w/at-4/results","Restricted y domain length with default settings",quant="THETA")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/orlanski-w/warm-bland-2/results","Restricted y domain length with default settings",quant="THETA")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/randtopo/tcline-at-glib-ct2/results","Restricted y domain length with default settings",quant="THETA")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/orlanski-w/at-4/results","Restricted y domain length with default settings",quant="THETA")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/orlanski-w/at/results","Restricted y domain length with default settings",quant="THETA")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/squish-polyna/under/results","Restricted y domain length with default settings",quant="SALT")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG/results","Restricted y domain length with default settings")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-steep/results","Restricted y domain length with default settings")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-2O-rand/results","Restricted y domain length with default settings")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-2O/results","Restricted y domain length with default settings")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/tcline/at/results","Restricted y domain length with default settings")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/tcline/under/results","Restricted y domain length with default settings")

# fig,axises = plt.subplots(2,2)
# timeSeriesDashboard("/home/garrett/Projects/MITgcm_ISC/experiments/squish-polyna/above/results","under",fig,axises)
# plt.show()
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/randtopo/crashtest1/results","Restricted y domain length with default settings")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/randtopo/crashtest-orlanskiw/results","Restricted y domain length with default settings")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/randtopo/tcline-at-glib/results","Restricted y domain length with default settings")
bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG/results","Restricted y domain length with default settings",quant="SALT")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-2O/results","Restricted y domain length with default settings",quant="THETA")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-2O-rand/results","Restricted y domain length with default settings",quant="THETA")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-steep/results","Restricted y domain length with default settings",quant="THETA")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-2O/results","Restricted y domain length with default settings",quant="THETA")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-2O-rand/results","Restricted y domain length with default settings",quant="THETA")

#surfaceAnim("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG/results","Restricted y domain length with default settings",quant="THETA")
#surfaceAnim("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-2O/results","Restricted y domain length with default settings",quant="THETA")
#surfaceAnim("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-2O-rand/results","Restricted y domain length with default settings",quant="THETA")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/reference/PIG-steep/results","Restricted y domain length with default settings",quant="THETA")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/tcline/under/results","Restricted y domain length with default settings",quant="THETA")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/tcline/at/results","Restricted y domain length with default settings",quant="THETA")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/tcline/above/results","Restricted y domain length with default settings",quant="THETA")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/orlanski-w/warm-bland/results","Restricted y domain length with default settings",quant="THETA")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/randtopo/tcline-at-glib-ct2/results","Restricted y domain length with default settings",quant="THETA")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/squish_polyna/at/results","Restricted y domain length with default settings")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/squish_polyna/above/results","Restricted y domain length with default settings")
#timeSeriesDashboard("/home/garrett/Projects/MITgcm_ISC/experiments/tcline/under/results","Lower thermocline 4km resolution")
