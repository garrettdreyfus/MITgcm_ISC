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
import pickle
from jmd95 import dens
from tabulate import tabulate

#matplotlib.use("TkAgg")

def timeSeries(fname):

    fnameparts = fname.split("/")
    slug = fnameparts[-3]+"-"+fnameparts[-2]+".pickle"

    if os.path.isfile("data/modelTimeSeries/"+slug) and True:
        with open("data/modelTimeSeries/"+slug,"rb") as f:
            output = pickle.load(f)
        return output


    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,prefix=["THETA","SALT","momKE","SHIfwFlx","VVEL","UVEL","WVEL"],ignore_unknown_vars=True,iters=times,extra_variables = extra_variables)
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

    icemaskm = np.logical_and(icedraft>0,icedraft<h)
    for k in range(shflx.shape[0]):
        shiwflxs.append(np.mean(shflx[k][icemaskm.T])*(60*60*24*365)*(1/920.0))


    avgbts = barotropic_streamfunction_max(fname)
    print(fname)
    print(shiwflxs)
    nanmask = ~np.isnan(shiwflxs)
    #np.sum(shflx*mask,axis=(1,2))*(60*60*24*365)*(1/920.0)*(1/np.sum(mask,axis=(1,2)))
    bottomtemps = thetas    # print(bottomtemps.shape)
    icem = icemask(fname,ds)
    print(THETA.shape,icem.shape)
    print(icem.shape,icemaskm.shape)
    THETA = ds.THETA.values
    SALT = ds.SALT.values
    UVEL = ds.UVEL.values
    VVEL = ds.VVEL.values
    WVEL = ds.WVEL.values
    VEL = np.sqrt(UVEL**2 + VVEL**2 + WVEL**2)
    #plt.imshow(icemaskm)
    #plt.show()
    icem[:,~icemaskm.T] = 0
    icesurfacetemps = []
    icesurfacevels = []
    icesurfacesalts = []
    meltapprox=[]

    for k in range(THETA.shape[0]):
        icesurfacetemps.append(np.nansum(THETA[k][icem])/np.sum(icem))
        icesurfacevels.append(np.nansum(VEL[k][icem])/np.sum(icem))
        icesurfacesalts.append(np.nansum(SALT[k][icem])/np.sum(icem))
        meltapprox.append(np.nansum((THETA[k]*VEL[k])[icem])/np.sum(icem))

    gprimes = fullOnGPrime(ds,fname)
    print("GPRIMES",gprimes)
    output = {"ts":np.asarray(ts)[nanmask],"theta":np.asarray(thetas)[nanmask],"salt":np.asarray(salts)[nanmask],\
        "kes":np.asarray(kes)[nanmask],"avgbts":np.asarray(avgbts)[nanmask],\
        "shiflx":np.asarray(shiwflxs)[nanmask],"bottemp":np.asarray(bottomtemps)[nanmask],"icesurfacetemp":np.asarray(icesurfacetemps)[nanmask],\
        "icesurfacesalt":np.asarray(icesurfacesalts)[nanmask],\
        "icesurfacevel":np.asarray(icesurfacevels)[nanmask],"meltapprox":np.asarray(meltapprox)[nanmask], "incavity":np.asarray(incavity)[nanmask],"gprime":np.asarray(gprimes)}
    with open("data/modelTimeSeries/"+slug,"wb") as f:
        pickle.dump(output,f)
    return output

def matVarsFile(fname):
    one_up = os.path.dirname(fname)
    return one_up+"/input/metaparameters.mat"

def grabMatVars(fname,val_tup):
    variables = scipy.io.loadmat(matVarsFile(fname),variable_names=val_tup,struct_as_record=False)
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
    result = quad(f_interp,depth,min(depth+25,0), points = zz[::-1])[0]
    result = ((result/min(25,abs(depth))))
    return result

def maxTemp(depth,fname):
    print(depth)
    variables = grabMatVars(fname,('tNorth','tEast','zz',''))
    tEast = np.asarray(variables["tEast"])#[0]+1.8
    zz = np.asarray(variables["zz"])[0]
    tEast = tEast[int(tEast.shape[0]-1)]+1.8
    return zz[np.argmax(tEast)]



def slope(zglib,fname):
    variables = grabMatVars(fname,('Zcdw_pt_shelf','icedraft','tEast','zz','yy',"xx","Yicefront"))
    zpyc = np.asarray(variables["Zcdw_pt_shelf"])[0][0]
    yy = np.asarray(variables["yy"])[0]
    xx = np.asarray(variables["xx"])[0]
    XX,YY = np.meshgrid(xx,yy)
    y = np.asarray(variables["Yicefront"])[0][0]
    icedraft = np.asarray(variables["icedraft"])
    zgl = np.nanmin(icedraft)-np.nanmax(icedraft[icedraft!=0])
    print("slopes",zgl/y)
    return (abs(zgl)-200)/y , abs(zglib-zgl)/y

def volume(fname):
    variables = grabMatVars(fname,('icedraft',"h","YY","xx","yy","Yicefront"))
    icedraft = np.asarray(variables["icedraft"])
    h = np.asarray(variables["h"])
    YY = np.asarray(variables["YY"])
    yy = np.asarray(variables["yy"])[0]
    xx = np.asarray(variables["xx"])[0]
    XX,YY = np.meshgrid(xx,yy)
    diff = np.abs(h)-np.abs(icedraft)
    zgl = np.nanmin(icedraft)
    y = np.asarray(variables["Yicefront"])[0][0] - np.nanmin(YY.T[diff>0])
    grad = np.gradient(icedraft)
    grad[0] = (grad[0]/np.gradient(xx)[10])**2
    grad[1] = (grad[1]/np.gradient(yy)[10])**2
    grad = np.sqrt(grad[0] + grad[1])
    return np.nanmedian(grad[np.logical_and(icedraft!=0,diff!=0)])#np.mean(diff[np.logical_and(icedraft!=0,diff!=0)]) #+ abs(zglib-zgl)/y

def steadyStateAverage(fname,xval,fig,axises,color="blue",marker="o",title=""):
    ((ax1,ax2),(ax3,ax4)) = axises 
    data = timeSeries(fname)
    glib = GLIBfromFile(matVarsFile(fname))
    aisf = aisfdepth(matVarsFile(fname))
    #xval =0 
    for k in data.keys():
        if k != "ts":
            try:
                data[k] = np.nanmean(data[k][data["ts"]>2])
            except:
                data[k]=np.nan
    if "gprime" not in data.keys():
        data["gprime"] = np.nan
    variables = grabMatVars(fname,("Hshelf","Xeast","Xwest","randtopog_height","Yicefront","Zcdw_pt_South","tEast","sEast","icedraft","zz"))
    y = np.asarray(variables["Yicefront"])[0][0]
    icedraft = np.asarray(variables["icedraft"])

    lightdens = dens(np.asarray(data["icesurfacesalt"]),np.asarray(data["icesurfacetemp"]),0)
    zz = np.asarray(variables["zz"])[0]
    tNorth = np.asarray(variables["tEast"])[-1,:]
    sNorth = np.asarray(variables["sEast"])[-1,:]
    zgl = np.nanmin(icedraft)

    aisfxval = intTemp(aisf,fname)
    glibxval = intTemp(glib,fname)
    ices,beds = slope(aisf,fname)
    ices = volume(fname)

    #glibxval = iceFaceTemp(glib,fname)

    shelf_depth = variables["Hshelf"]
    shelf_width = variables["Xeast"]-variables["Xwest"]
    randtopog_height = variables["randtopog_height"]
    max_height = variables["Zcdw_pt_South"][0][0]
    ax1.scatter(shelf_width,lightdens,c=color,marker=marker)
    ax1.set_xlabel("Cavity Width AKA Lx")
    ax1.set_ylabel("sigma0 at ice interface")
    ## salt plot
    print(data)
    ax2.scatter(shelf_width,data["gprime"],c=color,marker=marker)
    ax2.set_xlabel("")
    ax2.set_ylabel("")
    deltaH = ((((max_height-75)/2.0+75) - glib))
    tcline_height = (max_height-75)/2.0+75
    ax3.scatter(aisfxval*np.abs(aisfxval)*ices,-data["shiflx"],c=color,marker=marker)
    ax3.set_xlabel("(off shore t at aisf) * (thermocline height above aisf)")
    ax3.set_ylabel("Melt Rate m/yr")
    
    #ax4.scatter(glib,-data["shiflx"],c=color,marker=marker)
    zpyci = np.argmin(np.abs(np.abs(zz)-abs(tcline_height)))
    localdens = dens(sNorth,tNorth,abs(tcline_height))
    #gprime = (localdens[zpyci-10] - localdens[zpyci+10])/(zz[zpyci-10] - zz[zpyci+10])
    gprime = (np.mean(localdens[zpyci:])-np.mean(localdens[:zpyci]))
    print(lightdens)
    #print(gprime)
    ax4.scatter((glibxval)*ices*(deltaH)*gprime,-data["shiflx"],c=color,marker=marker)
    #ax4.scatter((glibxval)*tcline_height*ices,-data["shiflx"],c=color,marker=marker)
    ax4.set_xlabel("(off shore t at hub) * (thermocline height above hub)")
    ax4.set_ylabel("Melt Rate m/yr")

def steadyStateAverageSimple(fname,xval,fig,ax1,title="",color="blue",marker="o"):
    data = timeSeries(fname)
    glib = GLIBfromFile(matVarsFile(fname))
    aisf = aisfdepth(matVarsFile(fname))
    #xval =0 

    for k in data.keys():
        if k != "ts":
            try:
                data[k] = np.nanmean(data[k][data["ts"]>2])
            except:
                data[k]=np.nan
    if "gprime" not in data.keys():
        data["gprime"] = np.nan

    variables = grabMatVars(fname,("Hshelf","Xeast","Xwest","randtopog_height","Yicefront","Zcdw_pt_South","tEast","sEast","icedraft","zz","h"))
    y = np.asarray(variables["Yicefront"])[0][0]

    icedraft = np.asarray(variables["icedraft"])
    h = np.asarray(variables["h"])

    icemask = np.logical_and(icedraft!=0,np.abs(np.abs(icedraft)-np.abs(h))>5)

    lightdens = dens(np.asarray(data["icesurfacesalt"]),np.asarray(data["icesurfacetemp"]),0)
    zz = np.asarray(variables["zz"])[0]
    tNorth = np.asarray(variables["tEast"])[-1,:]
    sNorth = np.asarray(variables["sEast"])[-1,:]
    rho0 = dens(sNorth,tNorth,0)
    zgl = np.nanmin(icedraft)

    aisfxval = intTemp(aisf,fname)
    glibxval = intTemp(glib,fname)
    ices,beds = slope(aisf,fname)
    ices = volume(fname)

    ax1.set_xlabel("Thermocline height above HUB")
    ax1.set_ylabel("Melt Rate m/yr")
    plt.title(title)
    #plt.xlim(-350,0)
    #plt.ylim(0,25)
    shelf_depth = variables["Hshelf"]
    shelf_width = variables["Xeast"]-variables["Xwest"]
    randtopog_height = variables["randtopog_height"]
    max_height = variables["Zcdw_pt_South"][0][0]
    tcline_height = (max_height-75)/2.0+75

    zpyci = np.argmin(np.abs(np.abs(zz)-abs(tcline_height)))
    localdens = dens(sNorth,tNorth,abs(tcline_height))
    gprime = (np.mean(localdens[zpyci:])-np.mean(localdens[:zpyci]))/np.mean(localdens)

    deltaH = -(abs(tcline_height)- abs(glib))
    print(abs(deltaH))
    ax1.scatter((glibxval)*glibxval*ices,-data["shiflx"],c=color,marker=marker)
    #ax1.scatter(tcline_height-glib,-data["shiflx"],c=color,marker=marker)




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
def barotropic_streamfunction_graph(fname,description,times=np.array([]),res=1):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    shortname, fpath = outPath(fname) 
    moviewriter = FFMpegFileWriter(fps=1)
    fig,ax1 = plt.subplots()
    depth = depthFromdZ(ds)
    UFULL = ds.UVEL.values#
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
            frame = ax1.pcolormesh(-melt,cmap="jet",vmax=0.0012)
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
    inputname = "/".join(fname.split("/")[:-1])
    bmask = bottomMask(inputname,ds)
    icem = icemask(inputname,ds)
    shortname, fpath = outPath(fname) 
    moviewriter = FFMpegFileWriter(fps=1)
    fig,(ax1,ax2) = plt.subplots(1,2)
    theta = ds.THETA.values
    uvel = ds.UVEL.values
    vvel = ds.VVEL.values
    mask = np.full_like(uvel[0,0],1,dtype=bool)
    mask[::4,::4] = 0
    X,Y = np.meshgrid(range(uvel.shape[2]),range(uvel.shape[3]))
    print(fpath)
    with moviewriter.saving(fig, fpath+"-bottomvtop.mp4" , dpi=250):
        for k in tqdm(range(int(ds.UVEL.values.shape[0]))):
            thetak = theta[k]
            thetabot = thetak *bmask
            thetabot[~bmask] = np.nan
            thetabot = np.nanmean(thetabot,axis=0)
            ubot =  np.nanmean(uvel[k] * bmask,axis=0)
            vbot =  np.nanmean(vvel[k] * bmask,axis=0)
            ubot[mask] = np.nan
            vbot[mask] = np.nan

            utop =  np.nanmean(uvel[k] * icem,axis=0)
            vtop =  np.nanmean(vvel[k] * icem,axis=0)
            utop[mask] = np.nan
            vtop[mask] = np.nan
            frame = ax1.pcolormesh(thetabot,vmin=-2,vmax=2,cmap="jet")
            ax1.quiver(ubot,vbot,scale=0.05)
            ax1.set_title("Bottom")
            ax2.set_title("Top")

            thetatop = thetak *icem
            thetatop[~icem] = np.nan
            thetatop = np.nanmean(thetatop,axis=0)
            frame = ax2.pcolormesh(thetatop,vmin=-2,vmax=2,cmap="jet")
            if k==1:
                cb = plt.colorbar(frame,ax=ax2,pad=0.05)


            ax2.quiver(utop,vtop,scale=0.05)
            #plt.show()
            moviewriter.grab_frame()
            #cb.remove()
            frame.remove()

def iceFaceMelt(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    bmask = bottomMask(fname,ds)
    icem = icemask(fname,ds)
    shortname, fpath = outPath(fname) 
    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
    theta = ds.THETA.values
    uvel = ds.UVEL.values
    vvel = ds.VVEL.values
    wvel = ds.WVEL.values
    mask = np.full_like(uvel[0,0],1,dtype=bool)
    mask[::4,::4] = 0
    shifwflx = ds.SHIfwFlx.values

    X,Y = np.meshgrid(range(uvel.shape[2]),range(uvel.shape[3]))
    for k in tqdm(range(int(ds.UVEL.values.shape[0]))):
        thetak = theta[k]
        thetabot = thetak *bmask
        thetabot[~bmask] = np.nan
        thetabot = np.nanmean(thetabot,axis=0)
        utop =  np.nanmean(uvel[k] * icem,axis=0)
        vtop =  np.nanmean(vvel[k] * icem,axis=0)
        wtop =  np.nanmean(wvel[k] * icem,axis=0)
        frame = ax1.pcolormesh(utop)
        if k==1:
            cb = plt.colorbar(frame,ax=ax1,pad=0.1)
        #ax1.quiver(ubot,vbot,scale=1.5)

        thetatop = thetak *icem
        thetatop[~icem] = np.nan
        thetatop = np.nanmean(thetatop,axis=0)+1.8
        frame = ax2.pcolormesh(thetatop)
        if k==1:
            cb = plt.colorbar(frame,ax=ax2,pad=0.1)

        ax3.pcolormesh((thetatop**1)*np.sqrt(utop**2+vtop**2+wtop**2))
        ax4.pcolormesh(-shifwflx[k])
        plt.show()
        
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
    return np.logical_and(bottom_dist > -20,bottom_dist<0)


def fullOnGPrime(ds,fname):
    THETA = ds.THETA.values
    SALT = ds.SALT.values
    vals = grabMatVars(fname,("icedraft"))
    icedraft = np.abs(np.asarray(vals["icedraft"]))
    idmask=icedraft>0
    gprimes = []
    print("GPRIME TIME")
    for t_index in tqdm(range(THETA.shape[0])):
        gprimes_t = []
        for y_index in range(THETA.shape[2]):
            for x_index in range(THETA.shape[3]):
                tcast = THETA[t_index,:,y_index,x_index]
                if (np.sum(~np.isnan(tcast)))>10 and idmask[x_index,y_index]:
                    scast = SALT[t_index,:,y_index,x_index]
                    t = tcast[~np.isnan(scast)]
                    s = scast[~np.isnan(scast)]
                    d = dens(s,t,500)
                    if np.sum(d-d[0]>0.03)>0:#and np.sum(t>0.5)>0:
                        mldi = np.where(d-d[0]>0.03)[0][0]
                        #cdwi = np.where(t>0)[0][0]
                        rho_1 = np.nanmean(d[:mldi])
                        rho_2 = np.nanmean(d[mldi:])
                        gprimes_t.append((rho_2-rho_1)/np.mean((rho_1,rho_2)))
        gprimes.append(np.nanmean(gprimes_t))
    return gprimes




    
    
def outPath(resultspath):
    nameparts = resultspath.split("/")
    shortname = nameparts[-3] + "|" + nameparts[-2]
    fpath = "/".join((nameparts[:-5]+["pics"] + [shortname]))
    return shortname, fpath

def crossSectionAnim(fname,description,quant="THETA",res=1,dim="meridional"):
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
        if dim == "zonal":
            zonal_average = ds.where(ds.hFacC==1).mean(dim="XC",skipna=True)
            ys = zonal_average.YC.values
        if dim == "meridional":
            zonal_average = ds.isel(YC=50)
            ys = zonal_average.XC.values
        #zonal_average = ds.isel(XC=90)
        zvals = zonal_average[quant].values
        #zvals[zvals==0]=np.nan
        m = np.nanmedian(zvals)
        s = np.nanstd(zvals)
        tmin, tmax = m-5*s,m+s*2
        shortname, fpath = outPath(fname) 
        #plt.hist(zvals[:,:,:].flatten())
        #plt.show()
        fig.suptitle(shortname)
        zs = zonal_average.Z.values
        length = zvals.shape[0]
    if quant=="DENS":
        ds = open_mdsdataset(fname,prefix=["SALT","THETA"],ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
        #zonal_average = ds.where(ds.hFacC==1).isel(XC=100)
        if dim == "zonal":
            zonal_average = ds.where(ds.hFacC==1).mean(dim="XC",skipna=True)
            ys = zonal_average.YC.values
        if dim == "meridional":
            zonal_average = ds.isel(YC=50)
            ys = zonal_average.XC.values
        shortname, fpath = outPath(fname) 
        fig.suptitle(shortname)
        zs = zonal_average.Z.values
        #tmin, tmax = np.nanmin(zonal_average[quant]), np.nanmax(zonal_average[quant])
        zvals = (zonal_average["SALT"].values,zonal_average["THETA"].values)
        length = zvals[0].shape[0]
    print("starting movie")
    with moviewriter.saving(fig, fpath+quant+"|"+dim+".mp4" , dpi=250):
        for k in tqdm(range(0,length,res)):
            if quant == "DENS":
                frame = ax1.pcolormesh(ys,zs,gsw.sigma0(zvals[0][k,:,:],zvals[1][k,:,:]),cmap="jet",vmin=27.4,vmax=27.7)
            elif quant == "THETA":
                frame = ax1.imshow(zvals[k,:,:],cmap="jet",vmin=-0.5,vmax=1)
            elif quant == "SALT":
                frame = ax1.imshow(zvals[k,:,:],cmap="jet",vmin=34,vmax=34.8)
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
            nancount = icesurfacevelnp.nansum(np.isnan(znew),axis=0)
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
            #znew = np.multiply(zmask,X)gy
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


def generateRunsTable(fnames):
    table = []
    prettynames = {'shelf_depth':"Nominal depth of shelf", \
            'rng_seed':"Random bathymetry seed used", 'random_amplitude':"Amplitude of random bathymetry",\
            'cavity_depth':"Depth of cavity relative to depth of shelf", 'cavity_width': "Width of cavity",\
            'yicefront':"Northward extent of ice shelf"}#, 'tcline_atshelf_depth': "Depth of temperature maximum"}
    for fname in fnames:
        variables = grabMatVars(fname,("experiment_parameters"))
        fields = variables["experiment_parameters"][0][0].__dict__
        d = {}
        for l in fields.keys():
            if l in prettynames.keys():
                d[prettynames[l]] = fields[l][0][0]
        if d not in table:
            table.append(d)
    print(tabulate(table,headers="keys",tablefmt="latex",maxcolwidths=[8]*len(table[0].keys())))
    return table
        

#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-32/at125w100/results/","",dim="meridional",quant="DENS")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-32/at125w50/results/","",dim="meridional",quant="DENS")
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-32/at125w250/results/","",dim="meridional",quant="DENS")
#meltmap("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-32/at125w100/results/","")
#meltmap("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-32/at125w50/results/","")
#meltmap("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-32/at125w250/results/","")
#bottomAnim("/home/garrett/Projects/MITgcm_ISC/experiments/shelfzexp-GLIB-explore-101/at0d400/results/","")
# #meltmap("/home/garrett/Projects/MITgcm_ISC/experiments/shelfzexp-GLIB-explore-101/at0d400/results/","")
# fnames = [] 
# for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-32/*/", recursive = True):
#     fnames.append(f)
# for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/shelfzexp-GLIB-explore-101/*/", recursive = True):
#     fnames.append(f)
# for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/shelfzexp-highres-GLIB-explore-101/*/", recursive = True):
#     fnames.append(f)
# for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/icefront-GLIB-explore-32/*/", recursive = True):
#     fnames.append(f)
# for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/inverse-GLIB-explore-32/*/", recursive = True):
#     fnames.append(f)
# for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/slope200-GLIB-explore-18/*/", recursive = True):
#     fnames.append(f)
# for f in glob.glob("/home/garrett/Projects/MITgcm_ISC/experiments/slope375-GLIB-explore-18/*/", recursive = True):
#     fnames.append(f)
#
# generateRunsTable(fnames)
#
#


def folderMap(runsdict,save=True):
    fig,axises = plt.subplots(1,1)
    for k in runsdict.keys():
        print(k)
        for f in glob.glob(str("/home/garrett/Projects/MITgcm_ISC/experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                if key and key in f:
                    try:
                        steadyStateAverageSimple(f+"/results",0,fig,axises,color=runsdict[k]["color"][l],marker=runsdict[k]["marker"][l],title=runsdict[k]["description"][0])
                    except:
                        print("yeesh")
                elif not key:
                    try:
                        steadyStateAverageSimple(f+"/results",0,fig,axises,color=runsdict[k]["color"][l],marker=runsdict[k]["marker"][l],title=runsdict[k]["description"][0])
                    except:
                        print("yeesh")
            if save:
                plt.savefig("/home/garrett/Projects/HUB/paperfigures/"+k+".png")
        
runsdict = {\
        "widthexp-GLIB-explore-32":{"specialstring":['w50','w100','w250'], "marker":["o","o","o"] ,"color":["orange","purple","gray"],"description":["Different shelf widths" ]},\
        "shelfzexp-GLIB-explore-101":{"specialstring":['d200','d400','d500','d700','d800'], "marker":["x","x","x","x","x"] ,"color":["orange","purple","gray","green","red"],"description":["Different shelf depths"]},\
        "icefront-GLIB-explore-32":{"specialstring":['y100','y250','y280'], "marker":["^","^","^"] ,"color":["orange","purple","gray"],"description":["Different ice front distances"]},\
        "inverse-GLIB-explore-32":{"specialstring":['d0','d-200'], "marker":["D","D"] ,"color":["purple","orange"],"description":["Without a sill"]},\
        "slope200-GLIB-explore-18":{"specialstring":[False], "marker":["_"] ,"color":["orange"],"description":["Less steep slope"]},\
        "slope375-GLIB-explore-18":{"specialstring":[False], "marker":["_"] ,"color":["purple"],"description":["more steep slope"]}\
        }
folderMap(runsdict,save=False)

plt.show()


