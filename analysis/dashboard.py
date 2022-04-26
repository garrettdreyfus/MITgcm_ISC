import numpy as np
import glob
from xmitgcm import open_mdsdataset
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegFileWriter
from tqdm import tqdm



def timeSeriesDashboard(fname,label,fig,axises,times=np.array([])):
    ((ax1,ax2),(ax3,ax4)) = axises 
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))

    if times.any():
        ds = open_mdsdataset(fname,ignore_unknown_vars=True,iters=times,extra_variables = extra_variables)
    else:
        times=getIterNums(fname)
        ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    totalvolume  = (ds.XC*ds.YC*ds.Z*ds.hFacC).sum().values
    ## theta plot
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
    totalSHIfwFlx = (ds.SHIfwFlx*ds.XC*ds.YC*ds.Z*ds.hFacC).sum(axis=[1,2,3]).values
    ax4.plot(times,totalSHIfwFlx/totalvolume,label=label)
    ax4.set_xlabel("Time")
    ax4.set_ylabel("Fresh Water Flux")

def bottomMask(ds):
    bmask = np.full_like(ds.THETA[0,0,:,:],0,dtype=int)
    print(ds)
    for k in range(ds.Z.shape[0])[::-1]:
        nanmask = (bmask == 0)
        full_mask = np.full_like(bmask,k,dtype=int)
        full_mask[ds.maskC[k,:,:]==0] = 0
        bmask[nanmask] = full_mask[nanmask]
    return bmask
    

def crossSectionAnim(fname,description,times=np.array([]),quant="THETA"):
    fig,ax1 = plt.subplots()
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    print(times)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=list(times))
    ds[quant].values=ds[quant].values*ds.hFacC.values
    zonal_average = ds.where(ds.hFacC == 1).mean(dim="XC",skipna=True)
    #zonal_average = ds.isel(XC=32)
    moviewriter = FFMpegFileWriter(fps=1)
    tmin, tmax = np.nanmin(zonal_average[quant]), np.nanmax(zonal_average[quant])
    with moviewriter.saving(fig, 'myfile.mp4', dpi=250):
        for k in tqdm(range(zonal_average.THETA.shape[0])):
            frame = ax1.pcolormesh(zonal_average.YC.values,zonal_average.Z.values,zonal_average[quant][k,:,:],cmap="jet",vmin=tmin,vmax=tmax)
            cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            cb.remove()
            frame.remove()

def bottomAnim(fname,description,times=np.array([]),quant="THETA"):
    fig,ax1 = plt.subplots()
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    bmask = bottomMask(ds)

    #zonal_average = ds.isel(XC=32)
    moviewriter = FFMpegFileWriter(fps=1)
    tmin, tmax = np.nanmin(ds[quant]), np.nanmax(ds[quant])
    with moviewriter.saving(fig, 'myfile.mp4', dpi=250):
        print("writing movie")
        for k in tqdm([0]+list(range(ds[quant].shape[0]))+[-1]):
            frame = ax1.pcolormesh(ds.XC.values,ds.YC.values,ds[quant].values[k][bmask],cmap="jet",vmin=tmin,vmax=tmax)
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
    print(iters,saltiters)
    return np.unique(np.asarray(np.intersect1d(iters,saltiters)[:]))

#timeSeriesDashboard("/home/garrett/Projects/MITgcm_ISC/experiments/squish/test/results","Restricted y domain length with default settings",times = np.asarray(range(1,9))*420480)
# fig,axises = plt.subplots(2,2)
# timeSeriesDashboard("/home/garrett/Projects/MITgcm_ISC/experiments/tcline/under/results","under",fig,axises)
# timeSeriesDashboard("/home/garrett/Projects/MITgcm_ISC/experiments/tcline/above/results","above",fig,axises)
# timeSeriesDashboard("/home/garrett/Projects/MITgcm_ISC/experiments/tcline/at/results","at",fig,axises)
# plt.show()
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/squish-polyna/above/results","Restricted y domain length with default settings")
bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/squish-polyna/above/results","Restricted y domain length with default settings")
#timeSeriesDashboard("/home/garrett/Projects/MITgcm_ISC/experiments/tcline/under/results","Lower thermocline 4km resolution")
