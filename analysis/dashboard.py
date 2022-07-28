import numpy as np
import glob
from xmitgcm import open_mdsdataset
import numpy.matlib 
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegFileWriter
from tqdm import tqdm
import os
from pathlib import Path
import re
import scipy
from scipy.integrate import quad
from matlabglib import GLIBfromFile

def timeSeries(fname):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,prefix=["THETA","SALT","momKE","SHIfwFlx"],ignore_unknown_vars=True,iters=times,extra_variables = extra_variables)
    totalvolume  = 1#(ds.XC.diff("x")*ds.YC.diff("y")*ds.Z.diff("z")*ds.hFacC).sum().values
    ## theta plot
    ts = ds.time.values*grabDeltaT(fname)/60.0/60.0/24.0/365.0
    tsnew = np.full_like(ts,0,dtype=float)
    tsnew[:] = ts
    ts = tsnew
    ts = ts/1000000000
    times=ts
    yice = grabMatVars(fname,("Yicefront"))["Yicefront"][0][0]
    cavity = ds.where(ds.YC<yice,drop=True)
    volume = (ds.hFacC * ds.drF * ds.rA).values
    cavvolume = (cavity.hFacC * cavity.drF * cavity.rA).values
    volumesum = np.sum(volume)
    cavvolumesum = np.sum(cavvolume)
    thetas=[]
    THETA = cavity.THETA.values
    for k in range(cavity.THETA.shape[0]):
        thetas.append((np.sum(THETA[k]*cavvolume))/cavvolumesum)
    salts = []
    SALT = cavity.SALT.values
    for k in range(cavity.THETA.shape[0]):
        salts.append((np.sum(SALT[k]*cavvolume))/cavvolumesum)
    kes = []
    momKE = cavity.momKE.values
    for k in range(cavity.THETA.shape[0]):
        kes.append((np.sum(momKE[k]*cavvolume))/cavvolumesum)
    mask = ~np.isnan(ds.SHIfwFlx.values)
    shflx = ds.SHIfwFlx.values
    shiwflxs = np.sum(shflx*mask,axis=(1,2))*(60*60*24*365)*(1/920.0)*(1/np.sum(mask,axis=(1,2)))

    bottomtemps = thetas    # print(bottomtemps.shape)
    surfacetemps = thetas
    return {"ts":np.asarray(ts),"theta":np.asarray(thetas),"salt":np.asarray(salts),"kes":np.asarray(kes),"shiflx":np.asarray(shiwflxs),"bottemp":np.asarray(bottomtemps),"surfacetemp":np.asarray(surfacetemps)}

def matVarsFile(fname):
    one_up = os.path.dirname(fname)
    return one_up+"/input/metaparameters.mat"

def grabMatVars(fname,val_tup):
    variables = scipy.io.loadmat(matVarsFile(fname),variable_names=val_tup)
    return variables

def intTemp(depth,fname):
    variables = grabMatVars(fname,('tNorth','tEast','zz'))
    tEast = np.asarray(variables["tEast"])#[0]+1.8
    #tEast = tEast[int(tEast.shape[0]*(5/6))]+1.8
    tEast = tEast[-1]+1.8
    #tNorth = np.asarray(variables["tNorth"])[0]+1.8
    zz = np.asarray(variables["zz"])[0]
    f_interp = lambda xx: np.interp(xx, zz[::-1], tEast[::-1])
    print(depth)
    result = quad(f_interp,depth,min(depth+200,0), points = zz[::-1])[0]
    print(result/(min(200,abs(depth))))
    return result/(min(200,abs(depth)))

#def theAndrewFactio

def steadyStateAverage(fname,xval,fig,axises,color="blue"):
    ((ax1,ax2,ax5,ax7),(ax3,ax4,ax6,ax8)) = axises 
    data = timeSeries(fname)
    xval = intTemp(GLIBfromFile(matVarsFile(fname)),fname)
    for k in data.keys():
        if k != "ts":
            data[k] = np.nanmean(data[k][data["ts"]>2.5])
    variables = grabMatVars(fname,("Hshelf","randtopog_height","Zcdw_pt_South"))
    tcline_height = variables["Zcdw_pt_South"]
    shelf_depth = variables["Hshelf"]
    randtopog_height = variables["randtopog_height"]
    tcline_height = variables["Zcdw_pt_South"]
    ax1.scatter(xval,data["theta"],c=color)
    ax1.set_xlabel("Height above HUB")
    ax1.set_ylabel("Potential Temperature")
    ## salt plot
    ax2.scatter(xval,data["salt"],c=color)
    ax2.set_xlabel("Height above HUB")
    ax2.set_ylabel("Salinity")

    ## kinetic energy plot
    ax3.scatter(xval,data["kes"],c=color)
    ax3.set_xlabel("Height above HUB")
    ax3.set_ylabel("Kinetic Energy")

    ax4.scatter(xval,-data["shiflx"],c=color)
    ax4.set_xlabel("Height above HUB")
    ax4.set_ylabel("Melt Rate m/yr")
    
    ## kinetic energy plot
    #totalSHIfwFlx = (ds.SHIfwFlx*ds.XC*ds.YC*ds.Z*ds.hFacC).sum(axis=[1,2,3]).values
    bottomtemps = []
    surfacetemps = []

    ax5.scatter(xval,data["bottemp"],c=color)
    #ax6.scatter(xval,data["surfacetemp"],c=color)
    ax5.set_xlabel("Height above HUB")
    ax5.set_ylabel("Bottom Potential Temperature")
    #ax6.set_xlabel("Height above HUB")
    #ax6.set_ylabel("Surface Potential Temperature")
    ax7.scatter(shelf_depth,randtopog_height,c=color)
    ax7.set_xlabel("Depth of shelf before random perturbation")
    ax7.set_ylabel("Random perturbation amplitude")
    ax8.scatter(tcline_height,-data["shiflx"],c=color)
    ax8.set_xlabel("Depth of thermocline")
    ax8.set_ylabel("Melt Rate m/yr")


def timeSeriesDashboard(fname,label,fig,axises,times=np.array([])):
    ((ax1,ax2,ax5),(ax3,ax4,ax6)) = axises 
    data = timeSeries(fname)
    ax1.plot(data["ts"][data["ts"]>2.5],data["theta"][data["ts"]>2.5],label=label)
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Potential Temperature")
    ax1.legend()
    ## salt plot
    color = [l for l in ax1.lines if l._label == label][0]._color
    ax2.plot(data["ts"][data["ts"]>2.5],data["salt"][data["ts"]>2.5],label=label,c=color)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Salinity")

    ## kinetic energy plot
    kes = []
    ax3.plot(data["ts"][data["ts"]>2.5],data["kes"][data["ts"]>2.5],label=label,c=color)
    ax3.set_xlabel("Time")
    ax3.set_ylabel("Kinetic Energy")

    ax4.plot(data["ts"][data["ts"]>2.5],data["shiflx"][data["ts"]>2.5])
    ax4.set_xlabel("Time")
    ax4.set_ylabel("Melt Rate m/yr")
    
    ## kinetic energy plot
    #totalSHIfwFlx = (ds.SHIfwFlx*ds.XC*ds.YC*ds.Z*ds.hFacC).sum(axis=[1,2,3]).values
    bottomtemps = []
    surfacetemps = []

    ax5.plot(data["ts"][data["ts"]>2.5],data["bottemp"][data["ts"]>2.5],label=label + " bottom",c=color)
    ax6.plot(data["ts"][data["ts"]>2.5],data["surfacetemp"][data["ts"]>2.5],label=label + " surface",c=color)
    ax5.set_xlabel("Time")
    ax5.set_ylabel("Bottom Potential Temperature")
    ax6.set_xlabel("Time")
    ax6.set_ylabel("Surface Potential Temperature")
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
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

def grabDeltaT(fname):
    p = Path(fname).parents[0]
    p = str(p)+"/input/data"
    with open(p,"rb") as file:
        for line in file:
            s = str(line)
            if "deltaT" in s:
                result = re.search('=(.*),', s)
                return float(result.group(1))



def barotropic_streamfunction(fname,description,times=np.array([]),res=5):
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
        for k in tqdm(range(0,ds.UVEL.values.shape[0],res)):
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

def crossSectionAnim(fname,description,times=np.array([]),quant="THETA",res=1):
    fig,ax1 = plt.subplots()
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,prefix=quant,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    #ds[quant].values=ds[quant].values*ds.hFacC.values

    #print(ds.hFacC)
    zonal_average = ds.where(ds.hFacC ==1).mean(dim="XC",skipna=True)
    #zonal_average = ds.mean(dim="XC",skipna=True)
    #zonal_average = ds.isel(XC=50)
    moviewriter = FFMpegFileWriter(fps=1)
    tmin, tmax = np.nanmin(zonal_average[quant]), np.nanmax(zonal_average[quant])
    shortname, fpath = outPath(fname) 
    fig.suptitle(shortname)
    ys = zonal_average.YC.values
    zs = zonal_average.Z.values
    zvals = zonal_average[quant].values
    with moviewriter.saving(fig, fpath+".mp4" , dpi=250):
        for k in tqdm(range(0,zonal_average[quant].shape[0],res)):
            frame = ax1.pcolormesh(ys,zs,zvals[k,:,:],cmap="jet",vmin=tmin,vmax=tmax)
            cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            cb.remove()
            frame.remove()

def bottomAnim(fname,description,times=np.array([]),quant="SALT",res=5):
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
        for k in tqdm([0]+list(range(quantvals.shape[0]))+[-1]):
            d = quantvals[k]
            X = np.full_like(d,np.nan,dtype=float)
            X[ds.hFacC.values != 0]= d[ds.hFacC.values != 0]
            znew = np.multiply(zmask,X)
            nancount = np.nansum(np.isnan(znew),axis=0)
            znew = np.nansum(znew,axis=0)
            znew[nancount==X.shape[0]] = np.nan
            frame = ax1.pcolormesh(ds.XC.values,ds.YC.values,znew,cmap="jet",vmin=34.4,vmax=34.65)
            ax1.contour(ds.XC.values,ds.YC.values,depth,colors="black",levels=20)
            cb = plt.colorbar(frame)
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
        for k in tqdm([0]+list(range(0,quantvals.shape[0]))+[-1]):
            #X = np.full_like(quantvals[k],np.nan,dtype=float)
            #X[ds.hFacC.values != 0]= quantvals[k][ds.hFacC.values != 0]
            #znew = np.multiply(zmask,X)
            #nancount = np.nansum(np.isnan(znew),axis=0)
            #znew = np.nansum(znew,axis=0)
            #znew[nancount==X.shape[0]] = np.nan
            znew = quantvals[k][0,:,:]
            frame = ax1.pcolormesh(ds.XC.values,ds.YC.values,znew,cmap="jet",vmin=30)
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

# fig,axises = plt.subplots(2,3)
# timeSeriesDashboard("/home/garrett/Projects/MITgcm_ISC/experiments/saltflux-explore/up/results","up",fig,axises)
#fig,axises = plt.subplots(2,3)
#for k in [-200, -125, -50, 0, 50, 125, 200]:
    #try:
        #timeSeriesDashboard("/home/garrett/Projects/MITgcm_ISC/experiments/GLIB-explore-18/at"+str(k)+"/results",str(k),fig,axises)
    #except:
        #print("nope",k)
#plt.show()

# bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/tclinedz250-GLIB-explore-18/at-200/results","Restricted y domain length with default settings","SALT")
# bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/tclinedz250-GLIB-explore-18/at0/results","Restricted y domain length with default settings","SALT")
# crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/tclinedz250-GLIB-explore-18/at-200/results","Restricted y domain length with default settings","SALT")
# crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/tclinedz250-GLIB-explore-18/at0/results","Restricted y domain length with default settings","SALT")

fig,axises = plt.subplots(2,4)
for k in [-200, -125, -50, 0, 50, 125, 200]:
    steadyStateAverage("/home/garrett/Projects/MITgcm_ISC/experiments/GLIB-explore-16/at"+str(k)+"/results",k,fig,axises,color="red")

for k in [-200, -125, -50, 0, 50, 125, 200]:
    try:
        steadyStateAverage("/home/garrett/Projects/MITgcm_ISC/experiments/GLIB-explore-18/at"+str(k)+"/results",k,fig,axises,color="green")
    except:
        print("nope",k)
for k in [-200, -125, -50, 0, 50, 125, 200]:
    try:
        steadyStateAverage("/home/garrett/Projects/MITgcm_ISC/experiments/tclinedz250-GLIB-explore-18/at"+str(k)+"/results",k,fig,axises,color="blue")
    except:
        print("nope",k)
for k in [-200, -125, -50, 0, 50, 125, 200]:
    try:
        steadyStateAverage("/home/garrett/Projects/MITgcm_ISC/experiments/halfw-GLIB-explore-18/at"+str(k)+"/results",k,fig,axises,color="orange")
    except:
        print("nope",k)

for k in [-200, -125, -50, 0, 50, 125, 200]:
    try:
        steadyStateAverage("/home/garrett/Projects/MITgcm_ISC/experiments/doublew-GLIB-explore-18/at"+str(k)+"/results",k,fig,axises,color="gray")
    except:
        print("nope",k)

# for k in [-200, -100, -50, -25, 0, 150]:
#     try:
#         steadyStateAverage("/run/media/garrett/037e02f0-d92c-4b5e-8415-f3f936191171/experiments/GLIB-explore/at"+str(k)+"/results",k,fig,axises,color="purple")
#     except:
#         print("nope",k)
plt.show()
