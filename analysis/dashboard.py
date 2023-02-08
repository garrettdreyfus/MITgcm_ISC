import numpy as np
import glob
from xmitgcm import open_mdsdataset
import numpy.matlib 
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.animation import FFMpegFileWriter
from tqdm import tqdm
import os
from pathlib import Path
import re
import scipy
from scipy.integrate import quad
from matlabglib import GLIBfromFile,aisfdepth
import gsw

#matplotlib.use("TkAgg")

def timeSeries(fname):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,prefix=["THETA","SALT","momKE","SHIfwFlx","VVEL"],ignore_unknown_vars=True,iters=times,extra_variables = extra_variables)
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
    incavity = []
    vvel = ds.VVEL.values
    frontmask = ds.hFacC[:,np.nanargmin(np.abs(ds.VVEL.YG-yice)),:]==0
    sliceindex=np.nanargmin(np.abs(ds.VVEL.YG-yice))
    for k in range(ds.THETA.shape[0]):
        vvelcross = vvel[k,:,sliceindex,:]
        vvelcross[frontmask]=np.nan
        incavity.append(np.nansum(vvelcross[vvelcross<0]))
    mask = ~np.isnan(ds.SHIfwFlx.values)
    shflx = ds.SHIfwFlx.values
    shiwflxs = []
    vals = grabMatVars(fname,("h","icedraft"))
    icedraft = np.abs(np.asarray(vals["icedraft"]))
    h = np.abs(np.asarray(vals["h"]))
    icemask = np.logical_and(icedraft>0,icedraft<h)
    for k in range(shflx.shape[0]):
        shiwflxs.append(np.mean(shflx[k][icemask.T])*(60*60*24*365)*(1/920.0))
    avgbts = barotropic_streamfunction_max(fname)
    print(fname)
    print(shiwflxs)
    nanmask = ~np.isnan(shiwflxs)
    #np.sum(shflx*mask,axis=(1,2))*(60*60*24*365)*(1/920.0)*(1/np.sum(mask,axis=(1,2)))
    bottomtemps = thetas    # print(bottomtemps.shape)
    surfacetemps = thetas
    return {"ts":np.asarray(ts)[nanmask],"theta":np.asarray(thetas)[nanmask],"salt":np.asarray(salts)[nanmask],"kes":np.asarray(kes)[nanmask],"avgbts":np.asarray(avgbts)[nanmask],\
        "shiflx":np.asarray(shiwflxs)[nanmask],"bottemp":np.asarray(bottomtemps)[nanmask],"surfacetemp":np.asarray(surfacetemps)[nanmask], "incavity":np.asarray(incavity)[nanmask]}

def matVarsFile(fname):
    one_up = os.path.dirname(fname)
    return one_up+"/input/metaparameters.mat"

def grabMatVars(fname,val_tup):
    variables = scipy.io.loadmat(matVarsFile(fname),variable_names=val_tup)
    return variables

def intTemp(depth,fname):
    print(depth)
    variables = grabMatVars(fname,('tNorth','tEast','zz',''))
    tEast = np.asarray(variables["tEast"])#[0]+1.8
    zz = np.asarray(variables["zz"])[0]
    tEast = tEast[int(tEast.shape[0]-1)]+1.8
    f_interp = lambda xx: np.interp(xx, zz[::-1], tEast[::-1])
    results = []
    ls = []
    result = quad(f_interp,depth,min(depth+200,0), points = zz[::-1])[0]
    print("square")
    result = ((result/min(200,abs(depth)))**2)
    return result

def AndrewsMetric(zglib,fname):
    variables = grabMatVars(fname,('Zcdw_pt_shelf','icedraft','tEast','zz'))
    zpyc = np.asarray(variables["Zcdw_pt_shelf"])[0][0]
    icedraft = np.asarray(variables["icedraft"])
    zgl = np.nanmin(icedraft)
    print(zglib,zpyc,zgl)
    return (zglib-zgl)*(zpyc-zgl)


def steadyStateAverage(fname,xval,fig,axises,color="blue",marker="o"):
    ((ax1,ax2),(ax3,ax4)) = axises 
    data = timeSeries(fname)
    glib = GLIBfromFile(matVarsFile(fname))
    aisf = aisfdepth(matVarsFile(fname))
    #xval = AndrewsMetric(glib,fname)
    #xval =0 
    for k in data.keys():
        if k != "ts":
            data[k] = np.nanmean(data[k][data["ts"]>2])
    variables = grabMatVars(fname,("Hshelf","randtopog_height","Zcdw_pt_South","icedraft"))
    icedraft = np.asarray(variables["icedraft"])
    zgl = np.nanmin(icedraft)
    print(zgl,glib)
    xval = intTemp(aisf,fname)
    tcline_height = variables["Zcdw_pt_South"]
    #xval = tcline_height
    #xval = tcline_height
    #xval = tcline_height - glib
    shelf_depth = variables["Hshelf"]
    randtopog_height = variables["randtopog_height"]
    tcline_height = variables["Zcdw_pt_South"]
    ax1.scatter(xval,data["theta"],c=color,marker=marker)
    ax1.set_xlabel("Height above HUB")
    ax1.set_ylabel("Potential Temperature")
    ## salt plot
    ax2.scatter(xval,data["salt"],c=color,marker=marker)
    ax2.set_xlabel("Height above HUB")
    ax2.set_ylabel("Salinity")

    ## kinetic energy plot
    ax3.scatter(xval,data["avgbts"],c=color,marker=marker)
    ax3.set_xlabel("Height above HUB")
    ax3.set_ylabel("Kinetic Energy")

    ax4.scatter(xval,-data["shiflx"],c=color,marker=marker)
    ax4.set_xlabel("Height above HUB")
    ax4.set_ylabel("Melt Rate m/yr")
    
    ## kinetic energy plot
    #totalSHIfwFlx = (ds.SHIfwFlx*ds.XC*ds.YC*ds.Z*ds.hFacC).sum(axis=[1,2,3]).values
    # bottomtemps = []
    # surfacetemps = []

    # ax5.scatter(xval,data["bottemp"],c=color,marker=marker)
    # #ax6.scatter(xval,data["surfacetemp"],c=color)
    # ax5.set_xlabel("Height above HUB")
    # ax5.set_ylabel("Bottom Potential Temperature")
    #ax6.set_xlabel("Height above HUB")
    #ax6.set_ylabel("Surface Potential Temperature")

def steadyStateAverageSimple(fname,xval,fig,ax1,color="blue",marker="o"):
    data = timeSeries(fname)
    glib = GLIBfromFile(matVarsFile(fname))
    #xval = AndrewsMetric(glib,fname)
    #xval =0 
    for k in data.keys():
        if k != "ts":
            data[k] = np.nanmean(data[k][data["ts"]>2])
    variables = grabMatVars(fname,("Hshelf","randtopog_height","Zcdw_pt_South","icedraft"))
    icedraft = np.asarray(variables["icedraft"])
    zgl = np.nanmin(icedraft)
    print(zgl,glib)

    tcline_height = variables["Zcdw_pt_South"]
    #xval = tcline_height
    xval = intTemp(glib,fname)
    #xval = tcline_height - glib
    shelf_depth = variables["Hshelf"]
    randtopog_height = variables["randtopog_height"]
    tcline_height = variables["Zcdw_pt_South"]

    ax1.scatter(xval,-data["shiflx"],c=color,marker=marker)
    ax1.set_xlabel("Depth of thermocline above HUB")
    ax1.set_ylabel("Melt Rate m/yr")
    
    ## kinetic energy plot
    #totalSHIfwFlx = (ds.SHIfwFlx*ds.XC*ds.YC*ds.Z*ds.hFacC).sum(axis=[1,2,3]).values
    # bottomtemps = []
    # surfacetemps = []

    # ax5.scatter(xval,data["bottemp"],c=color,marker=marker)
    # #ax6.scatter(xval,data["surfacetemp"],c=color)
    # ax5.set_xlabel("Height above HUB")
    # ax5.set_ylabel("Bottom Potential Temperature")
    #ax6.set_xlabel("Height above HUB")
    #ax6.set_ylabel("Surface Potential Temperature")



def timeSeriesDashboard(fname,label,fig,axises,times=np.array([])):
    ((ax1,ax2,ax5),(ax3,ax4,ax6)) = axises 
    data = timeSeries(fname)
    starttime = 2.5
    ax1.plot(data["ts"][data["ts"]>starttime],data["theta"][data["ts"]>starttime],label=label)
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Potential Temperature")
    ax1.legend()
    ## salt plot
    color = [l for l in ax1.lines if l._label == label][0]._color
    ax2.plot(data["ts"][data["ts"]>starttime],data["salt"][data["ts"]>starttime],label=label,c=color)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Salinity")

    ## kinetic energy plot
    kes = []
    ax3.plot(data["ts"][data["ts"]>starttime],data["kes"][data["ts"]>starttime],label=label,c=color)
    ax3.set_xlabel("Time")
    ax3.set_ylabel("Kinetic Energy")

    ax4.plot(data["ts"][data["ts"]>starttime],data["shiflx"][data["ts"]>starttime])
    ax4.set_xlabel("Time")
    ax4.set_ylabel("Melt Rate m/yr")
    
    ## kinetic energy plot
    #totalSHIfwFlx = (ds.SHIfwFlx*ds.XC*ds.YC*ds.Z*ds.hFacC).sum(axis=[1,2,3]).values
    bottomtemps = []
    surfacetemps = []

    ax5.plot(data["ts"][data["ts"]>starttime],data["bottemp"][data["ts"]>starttime],label=label + " bottom",c=color)
    ax6.plot(data["ts"][data["ts"]>starttime],data["surfacetemp"][data["ts"]>starttime],label=label + " surface",c=color)
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

def grabFandB(fname):
    p = Path(fname).parents[0]
    p = str(p)+"/input/data"
    with open(p,"rb") as file:
        for line in file:
            s = str(line)
            if "f0" in s:
                print(s)
                f = re.search('=(.*),', s)
            if "beta=" in s:
                print(s)
                B = re.search('=(.*),', s)
    return float(f.group(1)),float(B.group(1))

def barotropic_streamfunction_graph(fname,description,times=np.array([]),res=1):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    shortname, fpath = outPath(fname) 
    moviewriter = FFMpegFileWriter(fps=1)
    fig,ax1 = plt.subplots()
    depth = depthFromdZ(ds)
    UFULL = ds.UVEL.values
    VFULL = ds.VVEL.values
    cavity_mask = (ds.SHIfwFlx[0].values !=0)
    with moviewriter.saving(fig, fpath+"-bt.mp4" , dpi=250):
        for k in tqdm(range(int(ds.UVEL.values.shape[0]*3/4),ds.UVEL.values.shape[0],res)):
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
            m,s = np.nanmedian(bt[cavity_mask]),np.nanstd(bt[cavity_mask])
            vmin,vmax = m-3*s,m+3*s
            frame = ax1.pcolormesh(ds.XC.values,ds.YC.values,bt,vmin=vmin,vmax=vmax)
            ax1.contour(ds.XC.values,ds.YC.values,depth,levels=20,colors="black",vmin=vmin,vmax=vmax)
            #ax1.quiver(ds.XC.values,ds.YC.values,np.sum(U,axis=0),np.sum(V,axis=0))
            plt.title(str(k))
            cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            cb.remove()
            frame.remove()



def barotropic_streamfunction_max(fname,times=np.array([]),res=1):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    shortname, fpath = outPath(fname) 
    depth = depthFromdZ(ds)
    UFULL = ds.UVEL.values
    VFULL = ds.VVEL.values
    bts = []
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
        bts.append(np.nanmax(bt))
    return bts


def meltmap(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    print(ds)
    shortname, fpath = outPath(fname) 
    moviewriter = FFMpegFileWriter(fps=1)
    fig,ax1 = plt.subplots()
    mask = np.logical_and(ds.hFacC.values[0]==0,np.sum(ds.hFacC.values,axis=0)!=0)
    theta = ds.THETA.values
    with moviewriter.saving(fig, fpath+"-meltmap.mp4" , dpi=100):
        for k in tqdm(range(int(ds.UVEL.values.shape[0]))):
            thetak = theta[k]
            thetaavg = np.nanmean(thetak,axis=0)
            melt= ds.SHIfwFlx.values[k]
            melt[~mask]=np.nan
            thetaavg[~mask]=np.nan
            frame = ax1.pcolormesh(-melt,cmap="jet")
            if k==0:
                cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            #cb.remove()
            frame.remove()

def mixmap(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    print(ds)
    shortname, fpath = outPath(fname) 
    moviewriter = FFMpegFileWriter(fps=1)
    fig,ax1 = plt.subplots()
    mask = np.logical_and(ds.hFacC.values[0]==0,np.sum(ds.hFacC.values,axis=0)!=0)
    theta = ds.THETA.values
    with moviewriter.saving(fig, fpath+"-mixmap.mp4" , dpi=100):
        for k in tqdm(range(int(ds.UVEL.values.shape[0]))):
            thetak = theta[k]
            thetaavg = np.nanmean(thetak,axis=0)
            melt= ds.MXLDEPTH.values[k]
            melt[~mask]=np.nan
            thetaavg[~mask]=np.nan
            print(np.nanmin(melt),np.nanmax(melt))
            frame = ax1.pcolormesh(melt,cmap="jet")
            if k==0:
                cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            #cb.remove()
            frame.remove()



def bottomVtop(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    bmask = bottomMask(fname,ds)
    icem = icemask(fname,ds)
    shortname, fpath = outPath(fname) 
    moviewriter = FFMpegFileWriter(fps=1)
    fig,(ax1,ax2) = plt.subplots(1,2)
    theta = ds.THETA.values
    uvel = ds.UVEL.values
    vvel = ds.VVEL.values
    mask = np.full_like(uvel[0,0],1,dtype=bool)
    mask[::4,::4] = 0
    X,Y = np.meshgrid(range(uvel.shape[2]),range(uvel.shape[3]))
    with moviewriter.saving(fig, fpath+"-bottomvtop.mp4" , dpi=250):
        for k in tqdm(range(int(ds.UVEL.values.shape[0]))):
            thetak = theta[k]
            thetabot = thetak *bmask
            thetabot[~bmask] = np.nan
            thetabot = np.nanmean(thetabot,axis=0)
            ubot =  np.nanmean(uvel[k] * bmask,axis=0)
            vbot =  np.nanmean(vvel[k] * bmask,axis=0)
            ubot[mask] = 0
            vbot[mask] = 0
            frame = ax1.pcolormesh(thetabot,vmin=-2,vmax=2)
            if k==1:
                cb = plt.colorbar(frame,ax=ax1,pad=0.1)
            #ax1.quiver(ubot,vbot,scale=1.5)

            thetatop = thetak *icem
            thetatop[~icem] = np.nan
            thetatop = np.nanmean(thetatop,axis=0)
            frame = ax2.pcolormesh(thetatop,vmin=-2,vmax=2)
            if k==1:
                cb = plt.colorbar(frame,ax=ax2,pad=0.1)


            #ax2.quiver(utop,vtop,scale=0.05)
            #plt.show()
            moviewriter.grab_frame()
            #cb.remove()
            frame.remove()

def bottomMask(fname,ds,thresh=10):
    vals = grabMatVars(fname,("h","icedraft"))
    h = np.abs(np.asarray(vals["h"]))
    print(ds)
    Znew = np.zeros(ds.THETA.shape[1:])
    for l in range(Znew.shape[0]):
        Znew[l,:,:]=ds.Z.values[l]
    bottom_dist = Znew-(-h.T)
    return np.logical_and(bottom_dist < 50,bottom_dist>0)

def icemask(fname,ds,thresh=10):
    vals = grabMatVars(fname,("h","icedraft"))
    icedraft = np.abs(np.asarray(vals["icedraft"]))
    print(np.nanmin(icedraft))
    Znew = np.zeros(ds.THETA.shape[1:])
    for l in range(Znew.shape[0]):
        Znew[l,:,:]=ds.Z.values[l]
    bottom_dist = Znew-(-icedraft.T)
    return np.logical_and(bottom_dist > -50,bottom_dist<0)


    
    
def outPath(resultspath):
    nameparts = resultspath.split("/")
    shortname = nameparts[-3] + "|" + nameparts[-2]
    fpath = "/".join((nameparts[:-4]+["pics"] + [shortname]))
    return shortname, fpath

def crossSectionAnim(fname,description,quant="THETA",res=1):
    fig,ax1 = plt.subplots()
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    extra_variables = dict( KPPdiffS = dict(dims=["k","j","i"], attrs=dict(standard_name="KPPDIFFS", units="kg/m^3")))
    times=getIterNums(fname)
    #ds[quant].values=ds[quant].values*ds.hFacC.values

    #print(ds.hFacC)
    #zonal_average = ds.where(ds.hFacC==1).mean(dim="XC",skipna=True)
    #zonal_average = ds.mean(dim="XC",skipna=True)
    moviewriter = FFMpegFileWriter(fps=1)
    if quant!="DENS":
        ds = open_mdsdataset(fname,prefix=quant,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
        print(ds)
        zonal_average = ds.where(ds.hFacC==1).mean(dim="XC",skipna=True)
        #zonal_average = ds.isel(XC=90)
        zvals = zonal_average[quant].values
        #zvals[zvals==0]=np.nan
        m = np.nanmedian(zvals)
        s = np.nanstd(zvals)
        tmin, tmax = m-2*s,m+s*2
        shortname, fpath = outPath(fname) 
        #plt.hist(zvals[:,:,:].flatten())
        #plt.show()
        fig.suptitle(shortname)
        ys = zonal_average.YC.values
        zs = zonal_average.Z.values
    if quant=="DENS":
        ds = open_mdsdataset(fname,prefix=["SALT","THETA"],ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
        zonal_average = ds.where(ds.hFacC==1).isel(XC=100)
        shortname, fpath = outPath(fname) 
        fig.suptitle(shortname)
        ys = zonal_average.YC.values
        zs = zonal_average.Z.values
        #tmin, tmax = np.nanmin(zonal_average[quant]), np.nanmax(zonal_average[quant])
        zvals = (zonal_average["SALT"].values,zonal_average["THETA"].values)
    print("starting movie")
    with moviewriter.saving(fig, fpath+quant+".mp4" , dpi=250):
        for k in tqdm(range(0,zvals.shape[0],res)):
            if quant == "DENS":
                frame = ax1.imshow(ys,zs,gsw.sigma0(zvals[0][k,:,:],zvals[1][k,:,:]),cmap="jet",vmin=27.4,vmax=27.7)
            elif quant == "THETA":
                frame = ax1.imshow(zvals[k,:,:],cmap="jet",vmin=-0.5,vmax=1)
            elif quant == "SALT":
                frame = ax1.imshow(zvals[k,:,:],cmap="jet",vmin=tmin,vmax=tmax)
            else:
                frame = ax1.imshow(zvals[k,:,:],cmap="jet")
            cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            cb.remove()
            frame.remove()
    plt.close()

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
    quantvals = ds[quant].values#
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


def garrettsLittleProject(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    shortname, fpath = outPath(fname) 
    moviewriter = FFMpegFileWriter(fps=1)
    mask = np.logical_and(ds.hFacC.values[0]==0,np.sum(ds.hFacC.values,axis=0)!=0)
    f0,beta = grabFandB(fname)
    f = ds.YC.values*beta+f0

    vals = grabMatVars(fname,("h","icedraft"))
    h = np.abs(np.asarray(vals["h"]))-np.abs(np.asarray(vals["icedraft"]))
    dhdy = np.diff(h,axis=1,append=np.nan)

    Bt = beta - (f/h)*dhdy

    llhs = ((Bt*h)/f).T
    Bt=Bt.T

    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)

    ax1.set_title(r"$\hat{w}$ from vorticity")
    ax2.set_title(r"Freshwater flux")
    ax3.set_title(r"$\frac{(B_t * h)}{f}$")
    ax4.set_title(r"$W_{vel} output$")
    meltsum = []
    whatsum = []
    with moviewriter.saving(fig, fpath+"-glp.mp4" , dpi=250):
        for k in tqdm(range(ds.UVEL.values.shape[0])):
            if k>70 and k<120:
                V = ds.VVEL.values[k]
                T = ds.THETA.values[k]
                W = ds.WVEL.values[k]
                melt= ds.SHIfwFlx.values[k]
                Vhat = np.sum(V,axis=0)
                What = np.sum(W,axis=0)
                what = Vhat*llhs
                what[melt==0]=np.nan
                llhs[melt==0]=np.nan
                Bt[melt==0]=np.nan
                What[melt==0]=np.nan
                melt[melt==0]=np.nan
                s = np.nanstd(what)/4
                m = np.nanmedian(what)
                ax1.imshow(what,vmax=m+s,vmin=m-s)
                ax2.imshow(-melt)
                s = np.nanstd(llhs)/4
                m = np.nanmedian(llhs)
                ax3.imshow(llhs,vmax=m+s,vmin=m-s)
                s = np.nanstd(What)/2
                m = np.nanmedian(What)
                ax4.imshow(-What,vmax=m+s,vmin=m-s)
                #plt.show()
                #cb.remove()
                #frame.remove()

                whatsum.append(np.nansum(What))
                meltsum.append(np.nansum(melt))
                moviewriter.grab_frame()
                ax1.clear()
                ax2.clear()
                ax3.clear()
                ax4.clear()
    plt.show()
    plt.scatter(whatsum,meltsum)
    plt.show()

# crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/slope375-GLIB-explore-18/at-125/results","Restricted y domain length with default settings","RHOAnoma")
# crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/slope375-GLIB-explore-18/at0/results","Restricted y domain length with default settings","RHOAnoma")
# crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/slope375-GLIB-explore-18/at125/results","Restricted y domain length with default settings","RHOAnoma")

#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/doublew-GLIB-explore-18/at125/results","Restricted y domain length with default settings","THETA")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/GLIB-explore-21/at-400/results","Restricted y domain length with default settings","SALT")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/GLIB-explore-21/at-400/results","Restricted y domain length with default settings","SALT")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/doublew-GLIB-explore-18/at0/results","Restricted y domain length with default settings","SALT")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/slope200-GLIB-explore-18/at-125/results","Restricted y domain length with default settings","SALT")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/slope200-GLIB-explore-18/at0/results","Restricted y domain length with default settings","SALT")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/slope200-GLIB-explore-18/at125/results","Restricted y domain length with default settings","SALT")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/slope375-GLIB-explore-18/at-125/results","Restricted y domain length with default settings","SALT")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/slope375-GLIB-explore-18/at0/results","Restricted y domain length with default settings","SALT")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/shelfzexp-highres-GLIB-explore-101/at100d700/results","","THETA")
#meltmap("/home/garrett/Projects/MITgcm_ISC/experiments/shelfzexp-highres-GLIB-explore-101/at100d700/results","")
#exit()

#bottomVtop("/home/garrett/Projects/MITgcm_ISC/experiments/slope375-GLIB-explore-18/at-125/results","")
#bottomVtop("/home/garrett/Projects/MITgcm_ISC/experiments/slope375-GLIB-explore-18/at125/results","")
#bottomVtop("/home/garrett/Projects/MITgcm_ISC/experiments/slope200-GLIB-explore-18/at125/results","")
#mixmap("/home/garrett/Projects/MITgcm_ISC/experiments/slope375-GLIB-explore-18/at-125/results","")
#mixmap("/home/garrett/Projects/MITgcm_ISC/experiments/slope375-GLIB-explore-18/at125/results","")

#bottomVtop("/home/garrett/Projects/MITgcm_ISC/experiments/GLIB-explore-21/at-200/results","")
#meltmap("/home/garrett/Projects/MITgcm_ISC/experiments/icefront-GLIB-explore-32/at125y250/results","")
#meltmap("/home/garrett/Projects/MITgcm_ISC/experiments/icefront-GLIB-explore-32/at-300y250/results","")
#exit()

fig,axises = plt.subplots(2,2)
plt.ylim(0,19)
# for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/GLIB-explore-21/*/", recursive = True):
#     try:
#         steadyStateAverageSimple(f+"results",0,fig,axises,color="green")
#     except:
#         print("nope",f)
# plt.title("A reference case with a random bathymetry")
# plt.savefig("/home/garrett/Downloads/agutalk/1.png")

# for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/slope200-GLIB-explore-18/*/", recursive = True):
#     steadyStateAverageSimple(f+"results",0,fig,axises,color="blue")

# plt.title("Experiment with a smaller cavity slope and a different random bathymetry")
# plt.savefig("/home/garrett/Downloads/agutalk/2.png")

# for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/slope375-GLIB-explore-18/*/", recursive = True):
#     try:
#         steadyStateAverageSimple(f+"results",0,fig,axises,color="red")
#     except:
#         print("nope",f)
# plt.title("Experiment with a larger cavity slope")
# plt.savefig("/home/garrett/Downloads/agutalk/3.png")

# for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/GLIB-explore-21/*/", recursive = True):
#     try:
#         steadyStateAverageSimple(f+"results",0,fig,axises,color="green")
#     except:
#         print("nope",f)

# plt.title("Experiment with yet another random bathymetry")
# plt.savefig("/home/garrett/Downloads/agutalk/4.png")
# #steadyStateAverage("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-30/w50/"+"results",0,fig,axises,color="orange")
# #steadyStateAverage("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-30/w100/"+"results",0,fig,axises,color="purple")
# #steadyStateAverage("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-30/w250/"+"results",0,fig,axises,color="gray")

for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-32/*/", recursive = True):
    try:
        if "w50" in f:
            steadyStateAverage(f+"results",0,fig,axises,color="orange")
        if "w100" in f:
            steadyStateAverage(f+"results",0,fig,axises,color="purple")
        if "w250" in f:
            steadyStateAverage(f+"results",0,fig,axises,color="gray")
    except:
        print("nope",f)

# plt.title("Experiments exploring different widths")
# plt.savefig("/home/garrett/Downloads/agutalk/5.png")

for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/shelfzexp-GLIB-explore-101/*/", recursive = True):
    try:
        if "d200" in f:
            steadyStateAverage(f+"results",0,fig,axises,color="orange",marker="x")
        if "d400" in f:
            steadyStateAverage(f+"results",0,fig,axises,color="purple",marker="x")
        if "d800" in f:
            steadyStateAverage(f+"results",0,fig,axises,color="gray",marker="x")
    except:
        print("nope",f)

for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/shelfzexp-highres-GLIB-explore-101/*/", recursive = True):
    try:
        steadyStateAverage(f+"results",0,fig,axises,color="yellow",marker="s")
    except:
        print("nope",f)
# plt.title("Experiments exploring different baseline shelf depths with a different random bathymetry")
# plt.savefig("/home/garrett/Downloads/agutalk/6.png")



# for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/saltflux-explore/*/", recursive = True):
#     try:
#         steadyStateAverage(f+"results",0,fig,axises,color="pink",marker="|")
#     except:
#         print("nope",f)

for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/icefront-GLIB-explore-32/*/", recursive = True):
    try:
        if "y100" in f:
            steadyStateAverage(f+"results",0,fig,axises,color="red",marker="v")
        if "y250" in f:
            steadyStateAverage(f+"results",0,fig,axises,color="red",marker="^")
        if "y280" in f:
            steadyStateAverage(f+"results",0,fig,axises,color="red",marker="|")
    except:
        print("nope",f)




plt.show()


