import numpy as np
import glob
from xmitgcm import open_mdsdataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy.matlib 
import matplotlib.pyplot as plt
import matplotlib
import cmocean
from matplotlib.animation import FFMpegFileWriter
from sklearn.linear_model import LinearRegression
from tqdm import tqdm
import os
from pathlib import Path
import re
import scipy
from scipy.integrate import quad
from matlabglib import GLIBfromFile,aisfdepth
from scipy import stats
import gsw
import pickle
from jmd95 import dens
from tabulate import tabulate
from matplotlib.patches import Polygon

#matplotlib.use("TkAgg")

def timeSeries(fname,refresh=False):

    fnameparts = fname.split("/")
    slug = fnameparts[-3]+"-"+fnameparts[-2]+".pickle"

    if os.path.isfile("data/modelTimeSeries/"+slug) and not refresh:
        with open("data/modelTimeSeries/"+slug,"rb") as f:
            output = pickle.load(f)
        return output

    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,prefix=["THETA","SALT","momKE","SHIfwFlx","VVEL","UVEL","WVEL","RHOAnoma"],ignore_unknown_vars=True,iters=times,extra_variables = extra_variables)
    totalvolume  = 1#(ds.XC.diff("x")*ds.YC.diff("y")*ds.Z.diff("z")*ds.hFacC).sum().values
    ## theta plot
    ts = ds.time.values*grabDeltaT(fname)/60.0/60.0/24.0/365.0
    tsnew = np.full_like(ts,0,dtype=float)
    tsnew[:] = ts
    ts = tsnew
    ts = ts/1000000000
    times=ts
    yice = grabMatVars(fname,("Yicefront"))["Yicefront"][0][0]
    yshelf = grabMatVars(fname,("Yshelfbreak"))["Yshelfbreak"][0][0]
    ycoast = grabMatVars(fname,("Ycoast"))["Ycoast"][0][0]
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
    ht = vvel*(ds.THETA.values+1.8)*(ds.RHOAnoma.values+1000)
    frontmask = ds.hFacC[:,np.nanargmin(np.abs(ds.VVEL.YG-(yice-0))),:]==0
    sliceindex=np.nanargmin(np.abs(ds.VVEL.YG-(yice-0)))
    for k in range(ds.THETA.shape[0]):
        vvelcross = vvel[k,:,sliceindex,:]
        htcross = ht[k,:,sliceindex,:]*numpy.array(ds.drF.values)[:,None]
        htcross[frontmask]=np.nan
        vvelcross[frontmask]=np.nan
        incavity.append(np.nansum(htcross[vvelcross<0]))

    mask = ~np.isnan(ds.SHIfwFlx.values)
    shflx = ds.SHIfwFlx.values
    shiwflxs = []
    vals = grabMatVars(fname,("h","icedraft"))
    icedraft = np.abs(np.asarray(vals["icedraft"]))
    h = np.abs(np.asarray(vals["h"]))

    icemaskm = np.logical_and(icedraft>0,icedraft<h)
    for k in range(shflx.shape[0]):
        shiwflxs.append(np.mean(shflx[k][icemaskm.T])*(60*60*24*365)*(1/920.0))


    avgbts = [] #barotropic_streamfunction_max(fname)
    nanmask = ~np.isnan(shiwflxs)
    #np.sum(shflx*mask,axis=(1,2))*(60*60*24*365)*(1/920.0)*(1/np.sum(mask,axis=(1,2)))
    bottomtemps = thetas    # print(bottomtemps.shape)
    icem = icemask(fname,ds)
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

    gprimes,ssurf,scdw,tsurf,tcdw,froude = mixedLayerQuant(ds,fname)
    #gprimes,ssurf,scdw,tsurf,tcdw,froude = (np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)
    ###output = {"ts":np.asarray(ts)[nanmask],"theta":np.asarray(thetas)[nanmask],"salt":np.asarray(salts)[nanmask],\
        #"kes":np.asarray(kes)[nanmask],"avgbts":np.asarray(avgbts)[nanmask],\
        #"shiflx":np.asarray(shiwflxs)[nanmask],"bottemp":np.asarray(bottomtemps)[nanmask],"icesurfacetemp":np.asarray(icesurfacetemps)[nanmask],\
        #"icesurfacesalt":np.asarray(icesurfacesalts)[nanmask],\
        #"icesurfacevel":np.asarray(icesurfacevels)[nanmask],"meltapprox":np.asarray(meltapprox)[nanmask], "incavity":np.asarray(incavity)[nanmask],"gprime":np.asarray(gprimes)}
    output = {"ts":np.asarray(ts)[nanmask],"theta":np.asarray(thetas)[nanmask],"salt":np.asarray(salts)[nanmask],\
        "kes":np.asarray(kes)[nanmask],"avgbts":np.asarray([]),\
        "shiflx":np.asarray(shiwflxs)[nanmask],"bottemp":np.asarray(bottomtemps)[nanmask],"icesurfacetemp":np.asarray(icesurfacetemps)[nanmask],\
        "icesurfacesalt":np.asarray(icesurfacesalts)[nanmask],\
        "icesurfacevel":np.asarray(icesurfacevels)[nanmask],"meltapprox":np.asarray(meltapprox)[nanmask],\
        "incavity":np.asarray(incavity)[nanmask],"gprime":np.asarray(gprimes),\
        "ssurf":np.asarray(ssurf),"scdw":np.asarray(scdw),"tsurf":np.asarray(tsurf),"tcdw":np.asarray(tcdw),"froude":np.asarray(froude)\
    }
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
    #grad = np.sqrt(grad[0] + grad[1])
    grad = np.sqrt(np.sum(grad,axis=0))

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
                data[k] = np.nanmean(data[k][data["ts"]>5])
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
    #ax1.scatter(shelf_width,lightdens,c=color,marker=marker)
    #ax1.set_xlabel("Cavity Width AKA Lx")
    #ax1.set_ylabel("sigma0 at ice interface")
    ## salt plot
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
    glibi = np.argmin(np.abs(np.abs(zz)-abs(glib)))
    localdens = dens(sNorth,tNorth,abs(0))
    drhodz = np.diff(localdens)/np.diff(zz)
    zpyc = np.nanmean(zz[1:][drhodz>np.nanquantile(drhodz,0.90)])
    zpyci = np.nanargmin(np.abs(np.abs(zz)-abs(zpyc)))
    #gprime = (localdens[zpyci-10] - localdens[zpyci+10])/(zz[zpyci-10] - zz[zpyci+10])
    gprime = 9.8*(np.mean(localdens[zpyci:glibi])-np.mean(localdens[:zpyci]))/np.mean(localdens[:glibi])
    ax1.scatter(gprime*shelf_width,data["gprime"],c=color,marker=marker)
    #print(gprime)
    ax4.scatter((glibxval)*ices*(deltaH)*gprime,-data["shiflx"],c=color,marker=marker)
    #ax4.scatter((glibxval)*tcline_height*ices,-data["shiflx"],c=color,marker=marker)
    ax4.set_xlabel("(off shore t at hub) * (thermocline height above hub)")
    ax4.set_ylabel("Melt Rate m/yr")

def FStheory(fname,xval):
    data = timeSeries(fname)
    glib = GLIBfromFile(matVarsFile(fname))
    aisf = aisfdepth(matVarsFile(fname))
    #xval =0 

    for k in data.keys():
        if k != "ts":
            try:
                data[k] = np.nanmean(data[k][data["ts"]>2.5])
            except:
                data[k]=np.nan
    if "gprime" not in data.keys():
        data["gprime"] = np.nan

    variables = grabMatVars(fname,("Hshelf","Xeast","Xwest","randtopog_height","Yicefront","Zcdw_pt_South","Zcdw_pt_shelf","tEast","sEast","icedraft","zz","h","yy","xx"))
    icedraft = np.asarray(variables["icedraft"])
    h = np.asarray(variables["h"])
    #
    tNorth = np.asarray(variables["tEast"])[-1,:]
    sNorth = np.asarray(variables["sEast"])[-1,:]
    rho0 = dens(sNorth,tNorth,500)
    zgl = np.nanmin(icedraft)

    shelf_width = float(variables["Xeast"])-float(variables["Xwest"])
    aisfxval = intTemp(aisf,fname)
    glibxval = intTemp(glib,fname)
    ices,beds = slope(aisf,fname)
    ices = volume(fname)
    max_height = variables["Zcdw_pt_shelf"][0][0]
    tcline_height = (max_height-75)/2.0+75
    localdens = dens(sNorth,tNorth,abs(0))
    zz = np.asarray(variables["zz"])[0]
    #zpyci = np.argmax(np.abs(np.diff(localdens)/np.diff(zz)))
    gradd = np.abs(np.diff(localdens)/np.diff(zz))
    
    tcline_height=np.mean(zz[:-1][gradd>np.quantile(gradd,0.8)])#+75
    zpyci = np.argmin(np.abs(np.abs(zz)-abs(tcline_height)))
    #plt.plot(localdens,zz)
    #plt.axhline(y=tcline_height)
    #plt.axhline(y=(max_height-75)/2.0+75,color="red")
    #plt.show()
    #gprime_ext = 9.8*(np.mean(localdens[zpyci:])-np.mean(localdens[:min(zpyci*2,len(localdens)-1)]))/np.mean(localdens)
    #gprime_ext = 9.8*(np.mean(localdens[zpyci:min(zpyci*2,len(localdens)-1)])-np.mean(localdens[:zpyci]))/np.mean(localdens[:min(zpyci*2,len(localdens)-1)])
    d = localdens
    if np.sum(d-d[0]>0.03)>0:#and np.sum(t>0.5)>0:
        mldi = np.where(d-d[0]>0.03)[0][0]

        #cdwi = np.where(t>0)[0][0]
        rho_1 = np.nanmean(d[:mldi])
        rho_2 = np.nanmean(d[mldi:min(mldi*2,len(d)-1)])
        gprime_ext = 9.8*(rho_2-rho_1)/np.mean((rho_1,rho_2))
    #gprime_ext = 9.8*(np.mean(localdens[zpyci:min(zpyci+10,len(localdens)-1)])-np.mean(localdens[max(0,zpyci-10):zpyci]))/np.mean(localdens[:min(zpyci*2,len(localdens)-1)])
    rho_1i = np.logical_and(zz<zz[zpyci],zz>zz[zpyci]-zpyci)
    rho_2i = np.logical_and(zz<zz[zpyci]+zpyci,zz>zz[zpyci])
    print(rho_1i)
    print(rho_2i)
    gprime_ext = 9.8*(np.mean(localdens[rho_1i])-np.mean(localdens[rho_2i]))/np.mean(localdens[np.logical_or(rho_1i,rho_2i)])

    deltaH = -(abs(tcline_height)- abs(glib))
    if "reference" in fname and "at125" in fname:
        print(tcline_height)

    f = 1.3*10**-4
    rho0 = 1025
    rhoi = 910
    Cp = 4186
    If = 334000

    print(glibxval,deltaH,gprime_ext,f,ices)
    #return glibxval*deltaH*gprime_ext/f*ices,-data["shiflx"]
    #return glibxval*deltaH*(gprime_ext/625000)/f*ices,-data["shiflx"]
    #return (glibxval)*deltaH*(gprime_ext)/(f)*ices,-data["shiflx"]/(60*60*24*365)
    return (glibxval)*deltaH*(gprime_ext)/(f)*ices,-data["shiflx"]/(60*60*24*365)
    #return (glibxval)*glibxval*ices,-data["shiflx"]/(60*60*24*365)

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


def meltmapmovie(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
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
            frame = ax1.pcolormesh(melt,cmap="jet")
            if k==0:
                cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            #cb.remove()
            frame.remove()


    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))

    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    inputname = "/".join(fname.split("/")[:-1])
    bmask = bottomMask(inputname,ds)
    icem = icemask(inputname,ds)
    shortname, fpath = outPath(fname) 
    fig,ax1 = plt.subplots(1,1)
    theta = ds.THETA.values
    X = ds.XC.values
    Y = ds.YC.values
    thetatop = np.full_like(theta[1,1],0)
    count=0
    for k in tqdm(range(int(theta.shape[0]))):
        if times[k]>5:
            count+=1
            thetak = theta[k]
            thetatopc = thetak *icem
            thetatopc[~icem]=np.nan
            #thetatopc=thetak[0]
            if np.sum(~np.isnan(thetatopc)) !=0:
                thetatop += np.nanmean(thetatopc,axis=0)
    thetatop = thetatop/count
    mask = np.logical_and(ds.hFacC.values[0]==0,np.sum(ds.hFacC.values,axis=0)!=0)
    thetatop[~mask]=np.nan
    
    cax = ax1.pcolormesh(X/1000,Y/1000,thetatop,cmap=cmocean.cm.thermal,vmin=-0.8,vmax=0.2)
    cb = plt.colorbar(cax,ax=ax1,pad=0.05)

    plt.xlabel("x (km)",fontsize=18)
    plt.ylabel("y (km)",fontsize=18)
    plt.xlim(50,350)
    plt.ylim(0,150)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    #ax1.set_title("Width: 100",fontsize=18)
    plt.tight_layout
    plt.show()


def bottomVtop(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,prefix=["SALT","THETA","UVEL","VVEL","WVEL"],ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    inputname = "/".join(fname.split("/")[:-1])
    bmask = bottomMask(inputname,ds)
    icem = icemask(inputname,ds)
    shortname, fpath = outPath(fname) 
    fig,(ax1,ax2) = plt.subplots(1,2)
    theta = ds.THETA.values
    salt = ds.SALT.values
    uvel = ds.UVEL.values
    vvel = ds.WVEL.values
    X,Y = np.meshgrid(range(uvel.shape[3]),range(uvel.shape[2]))
    thetabot,saltbot,ubot,vbot = np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0)
    thetatop,salttop,utop,vtop = np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0)
    count=0
    for k in tqdm(range(int(ds.UVEL.values.shape[0]))):
        if times[k]>5:
            count+=1
            thetak = theta[k]
            thetabotc = thetak *bmask
            thetabotc[~bmask] = np.nan
            thetabot += np.nanmean(thetabotc,axis=0)

            saltk = salt[k]
            saltbotc = saltk *bmask
            saltbotc[~bmask] = np.nan
            #plt.imshow(saltbotc[:,:,90])
            #plt.show()
            #plt.imshow(ds.hFacC[:,:,90])
            #plt.show()
            saltbot += np.nanmean(saltbotc,axis=0)

            ubot += np.nanmean(uvel[k] * bmask,axis=0)
            vbot +=  np.nanmean(vvel[k] * bmask,axis=0)
            thetatopc = thetak *icem
            thetatopc[~bmask] = np.nan
            thetatop += np.nanmean(thetatopc,axis=0)
            utop +=  np.nanmean(uvel[k] * icem,axis=0)
            vtop +=  np.nanmean(vvel[k] * icem,axis=0)
    thetabot = thetabot/count
    thetatop = thetatop/count
    saltbot = saltbot/count
    cax = ax1.pcolormesh(X,Y,thetabot,vmin=-2,vmax=1,cmap=cmocean.cm.thermal)
    plt.colorbar(cax,ax=ax1)


    cax = ax2.pcolormesh(X,Y,vbot,vmin=-0.0001,vmax=0.0001,cmap="jet")
    #cax = ax2.pcolormesh(X,Y,saltbot,vmin=33.6,vmax=34.6,cmap=cmocean.cm.haline)
    plt.colorbar(cax,ax=ax2)
    
    ax1.set_title("Bottom temp")
    ax2.set_title("Bottom V")

    if k==1:
        cb = plt.colorbar(frame,ax=ax2,pad=0.05)

    plt.show()

def bottomspinup(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,prefix=["SALT","THETA","UVEL","VVEL","WVEL","RHOAnoma"],ignore_unknown_vars=True,extra_variables = extra_variables)
    inputname = "/".join(fname.split("/")[:-1])
    bmask = bottomMask(inputname,ds)
    icem = icemask(inputname,ds)
    shortname, fpath = outPath(fname) 
    theta = ds.THETA.values
    salt = ds.SALT.values
    uvel = ds.UVEL.values
    vvel = ds.WVEL.values
    X,Y = np.meshgrid(range(uvel.shape[3]),range(uvel.shape[2]))
    thetabot,saltbot,ubot,vbot = np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0)
    thetatop,salttop,utop,vtop = np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0)
    count=0
    for k in tqdm(range(10)):
        fig,(ax1,ax2) = plt.subplots(1,2)
        thetak = theta[k]
        thetak[~bmask]=np.nan
        thetabot = np.nanmean(thetak,axis=0)
        saltk = salt[k]
        saltk[~bmask]=np.nan
        saltbot = np.nanmean(saltk,axis=0)
        uvelk = uvel[k]
        uvelk[ds.hFacC==0]=np.nan
        ubot = np.nansum(uvelk,axis=0)
        vvelk = vvel[k]
        vvelk[ds.hFacC==0]=np.nan
        vbot = np.nansum(np.abs(vvelk),axis=0)
        cax = ax1.pcolormesh(X,Y,thetabot,vmin=-1,vmax=1,cmap=cmocean.cm.thermal)
        plt.colorbar(cax,ax=ax1)
        cax = ax2.pcolormesh(X,Y,dens(saltbot,thetabot,400),cmap=cmocean.cm.haline,vmin=1029.50,vmax=1029.75)
        #cax = ax2.pcolormesh(X,Y,vbot,cmap=cmocean.cm.haline,vmin=-0.005,vmax=0.005)
        plt.colorbar(cax,ax=ax2)
        fig.canvas.manager.full_screen_toggle() # toggle fullscreen mode
        plt.show()


    ax1.set_title("Bottom temp")
    ax2.set_title("Bottom V")

    if k==1:
        cb = plt.colorbar(frame,ax=ax2,pad=0.05)

    plt.show()


def someotherfunc():
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
        
    vals = grabMatVars(fname,("h","icedraft"))
    h = np.abs(np.asarray(vals["h"]))
    Znew = np.zeros(ds.THETA.shape[1:])
    for l in range(Znew.shape[0]):
        Znew[l,:,:]=ds.Z.values[l]
    bottom_dist = Znew-(-h.T)
    return np.logical_and(bottom_dist < 50,bottom_dist>0)


    
def outPath(resultspath):
    nameparts = resultspath.split("/")
    shortname = nameparts[-3] + "|" + nameparts[-2]
    fpath = "/".join((nameparts[:-5]+["pics"] + [shortname]))
    return shortname, fpath

def crossAndMelt(fname,name=""):
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(15,6))
    plt.subplots_adjust(wspace=0.45)
    crossSectionAverage(fname,"Reference",quant="THETA",dim="zonal",ax1=ax1,show=False)
    meltMapAverage(fname,"",ax1=ax2)

    plt.savefig("/home/garrett/Projects/HUB/paperfigures/crossAndMelts/"+name+".png")
    plt.close()



def generateRunsTable(runsdict):
    table = []
    prettynames = {'shelf_depth':"Nominal depth of shelf", \
            'rng_seed':"Random bathymetry seed used", 'random_amplitude':"Amplitude of random bathymetry",\
            'cavity_depth':"Depth of cavity relative to depth of shelf", 'cavity_width': "Width of cavity",\
            'yicefront':"Northward extent of ice shelf"}#, 'tcline_atshelf_depth': "Depth of temperature maximum"}
    for k in runsdict.keys():
        for f in glob.glob(str("/home/garrett/Projects/MITgcm_ISC/experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                fname = f+"/" #+"results/"
                variables = grabMatVars(fname,("experiment_parameters"))
                fields = variables["experiment_parameters"][0][0].__dict__
                conversion={"d0":"cavityd0","d-200":"cavityd-200","slope200":"cavityd200","slope375":"cavityd375"}
                d={}
                if not runsdict[k]["specialstring"][0]:
                    ss = k.split("-")[0]
                else:
                    ss = runsdict[k]["specialstring"][l]
                if ss in conversion.keys():
                    ss = conversion[ss]
                d["Experiment Name"] = ss
                for j in fields.keys():
                    if j in prettynames.keys():
                        d[prettynames[j]] = fields[j][0][0]
                if "Northward extent of ice shelf" not in d:
                    d["Northward extent of ice shelf"]=150
                if d not in table:
                    table.append(d)
    for k in table:print(k.keys())
    print(tabulate(table,headers="keys",tablefmt="latex",maxcolwidths=[9]*len(table[0].keys())))
    return table
        
def legendFunction(runsdict):
    conversion={"d0":"cavityd0","d-200":"cavityd-200","slope200":"cavityd200","slope375":"cavityd375"}
    for k in runsdict.keys():
        for l in range(len(runsdict[k]["specialstring"])):
            if not runsdict[k]["specialstring"][0]:
                ss = k.split("-")[0]
            else:
                ss = runsdict[k]["specialstring"][l]
            if ss in conversion.keys():
                ss = conversion[ss]
            plt.scatter(1,1,marker=runsdict[k]["marker"][l],color=runsdict[k]["color"][l],label=ss)
    plt.legend()
    plt.show()
            

def gprimeWidth(fname,xval,fig,ax1,title="",color="blue",marker="o"):
    #if "width" not in fname:# or 'at125' not in fname:
        #return 0
    data = timeSeries(fname)
    for k in data.keys():
        if k != "ts":
            try:
                data[k] = np.nanmean(data[k][data["ts"]>2])
            except:
                data[k]=np.nan
    if "gprime" not in data.keys():
        data["gprime"] = np.nan
    print(data["gprime"])
    variables = grabMatVars(fname,("Xeast","Xwest","Yicefront"))
    yice = float(variables["Yicefront"])
    shelf_width = float(variables["Xeast"])-float(variables["Xwest"])
    shiflx = -data["shiflx"]
    xval,shiflx = FStheory(fname,xval)
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))

    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    inputname = "/".join(fname.split("/")[:-1])
    icem = icemask(fname,ds)
    variables = grabMatVars(fname,("Hshelf","Xeast","Xwest","randtopog_height","Yicefront","Zcdw_pt_South","tEast","sEast","icedraft","zz","h","yy","xx"))
    icedraft = np.asarray(variables["icedraft"])
    h = np.asarray(variables["h"])
    #
    tNorth = np.asarray(variables["tEast"])[-1,:]
    sNorth = np.asarray(variables["sEast"])[-1,:]
    rho0 = dens(sNorth,tNorth,0)
    zgl = np.nanmin(icedraft)

    max_height = variables["Zcdw_pt_South"][0][0]
    tcline_height = (max_height-75)/2.0+75
    localdens = dens(sNorth,tNorth,abs(200))
    d = localdens
    Zfull = np.asarray(list(ds.Z))
    if np.sum(abs(tNorth-tNorth[0])>0.2)>0:#and np.sum(t>0.5)>0:
        t = tNorth[sNorth>0.1]
        Z = Zfull[sNorth>0.1]
        s = sNorth[sNorth>0.1]
        mldi = np.where(abs(tNorth-tNorth[0])>0.2)[0][0]
        d = dens(sNorth,tNorth,Z[mldi])
        #cdwi = np.where(t>0)[0][0]
        rho_1 = np.nanmean(d[:mldi])
        rho_2 = np.nanmean(d[mldi:min(mldi*2,len(d)-1)])
        #gprime_ext = np.nanmax(np.abs(np.diff(d)/np.diff(Z)))
        gprime_ext = 9.8*(rho_2-rho_1)/np.mean((rho_1,rho_2))

    #if np.sum(d-d[0]>0.03)>0:#and np.sum(t>0.5)>0:
        #mldi = np.where(d-d[0]>0.03)[0][0]

        #cdwi = np.where(t>0)[0][0]
        #rho_1 = np.nanmean(d[:mldi])
        #rho_2 = np.nanmean(d[mldi:min(mldi*5,len(d)-1)])
        #gprime_ext = 9.8*(rho_2-rho_1)/np.mean((rho_1,rho_2))
    zz = np.asarray(variables["zz"])[0]
    zpyci = np.argmin(np.abs(np.abs(zz)-abs(max_height)))
    #gprime_ext = 9.8*(np.mean(localdens[zpyci:min(zpyci*2,len(localdens)-1)])-np.mean(localdens[:zpyci]))/np.mean(localdens[:min(zpyci*2,len(localdens)-1)])
    f = 1.3*10**-4
    #ax1.scatter(gprime_ext*shelf_width,float(data["gprime"])*840,marker=marker,c=color)
    #ax1.scatter(gprime_ext*shelf_width,float(data["gprime"])*625000-525,marker=marker,c=color)
    #ax1.scatter(gprime_ext*shelf_width,float(data["gprime"]),marker=marker,c=color)
    ax1.scatter(gprime_ext*shelf_width,float(data["gprime"]),marker=marker,c=color)
    #plt.xlabel(r'$g^{\prime}_{ext}* W (m^2/s^2)$',fontsize=18)
    plt.xlabel(r'$g^{\prime}_{ext} (m/s^2)$',fontsize=18)
    plt.ylabel(r'$g^{\prime}_{in} (m/s^2)$',fontsize=18)
    #ax1.scatter(shelf_width,float(data["gprime"]),marker=marker,c=color)
    #ax1.scatter(shelf_width,gprime_ext,marker=marker,c=color)
    #ax1.plot(localdens,zz)
    #ax1.axhline(zz[zpyci])
    return gprime_ext

def gprimeAll(fname,xval,fig,ax1,title="",color="blue",marker="o"):
    data = timeSeries(fname)
    for k in data.keys():
        if k != "ts":
            try:
                data[k] = np.nanmean(data[k][data["ts"]>2])
            except:
                data[k]=np.nan
    if "gprime" not in data.keys():
        data["gprime"] = np.nan
    variables = grabMatVars(fname,("Xeast","Xwest","Yicefront"))
    yice = float(variables["Yicefront"])
    shelf_width = float(variables["Xeast"])-float(variables["Xwest"])
    shiflx = -data["shiflx"]
    xval,shiflx = FStheory(fname,xval)
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    extra_variables["KPPdiffT"] = dict(dims=["k","j","i"], attrs=dict(standard_name="KPP mld", units="m"))

    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    inputname = "/".join(fname.split("/")[:-1])
    icem = icemask(fname,ds)
    variables = grabMatVars(fname,("Hshelf","Xeast","Xwest","randtopog_height","Yicefront","Zcdw_pt_South","tEast","sEast","icedraft","zz","h","yy","xx"))
    icedraft = np.asarray(variables["icedraft"])
    h = np.asarray(variables["h"])
    #
    tNorth = np.asarray(variables["tEast"])[-1,:]
    sNorth = np.asarray(variables["sEast"])[-1,:]
    rho0 = dens(sNorth,tNorth,0)
    zgl = np.nanmin(icedraft)

    max_height = variables["Zcdw_pt_South"][0][0]
    tcline_height = (max_height-75)/2.0+75
    localdens = dens(sNorth,tNorth,abs(tcline_height))
    d = localdens
    Z = np.asarray(list(ds.Z))
    zz = np.asarray(variables["zz"])[0]
    if np.sum(abs(tNorth-tNorth[0])>0.2)>0:#and np.sum(t>0.5)>0:
        mldi = np.where(abs(tNorth-tNorth[0])>0.2)[0][0]
        #cdwi = np.where(t>0)[0][0]
        rho_1 = np.nanmean(d[:mldi])
        rho_2 = np.nanmean(d[mldi:min(mldi*2,len(d)-1)])
        #gprime_ext = np.nanmean(np.abs(np.diff(d)/np.diff(zz)))
        gprime_ext = 9.8*(rho_2-rho_1)/np.mean((rho_1,rho_2))

    zz = np.asarray(variables["zz"])[0]
    zpyci = np.argmin(np.abs(np.abs(zz)-abs(max_height)))
    gprime_ext = 9.8*(np.mean(localdens[zpyci:min(zpyci*2,len(localdens)-1)])-np.mean(localdens[:zpyci]))/np.mean(localdens[:min(zpyci*2,len(localdens)-1)])
    f = 1.3*10**-4
    #ax1.scatter(gprime_ext*shelf_width,float(data["gprime"])*840,marker=marker,c=color)
    #ax1.scatter(gprime_ext*shelf_width,float(data["gprime"])*625000-525,marker=marker,c=color)
    #ax1.scatter(gprime_ext*shelf_width,float(data["gprime"]),marker=marker,c=color)
    W0 = gprime_ext*shelf_width/float(data["gprime"])
    kpp = np.nanmean(ds.KPPdiffT.values,axis=0)
    kpp = np.nanmean(kpp,axis=0)
    kpp = np.nanmean(kpp[ds.YC.values<60000])
    small=np.log10(np.abs(ds.THETA_inst.values[1]))<-3
    smallmask=np.logical_and(small,ds.hFacC!=0)
    smallmask = np.sum(smallmask,axis=0)
    #plt.imshow(np.sum(smallmask,axis=0))
    #plt.title(fname)
    #plt.show()
    #XX,YY = np.meshgrid(ds.XC.values,ds.YC.values)
    #print("exact zeros",np.sum(np.log10(np.abs(ds.THETA.values[ds.hFacC!=0]))<-3))
    #print("exact zeros",np.sum(np.log10(np.abs(ds.THETA.values[ds.hFacC!=0]))<-3))
    #print("exact zeros",np.sum(np.log10(np.abs(ds.THETA.values[ds.hFacC!=0]))<-3))
    ax1.scatter(float(data["gprime"]),np.sum(smallmask[YY<60000]),marker=marker,c=color)
    #ax1.scatter(float(data["gprime"]),-float(data["shiflx"]),marker=marker,c=color)
    #plt.xlabel(r'$g^{\prime}_{in} (m/s^2)$',fontsize=18)
    #plt.ylabel(r'$m (m/yr)$',fontsize=18)
    #plt.xlabel(r'$g^{\prime}_{ext}* W (m^2/s^2)$',fontsize=18)
    plt.xlabel(r'$g^{\prime}_{ext} (m/s^2)$',fontsize=18)
    plt.ylabel(r'$g^{\prime}_{in} (m/s^2)$',fontsize=18)
    #ax1.scatter(shelf_width,float(data["gprime"]),marker=marker,c=color)
    #ax1.scatter(shelf_width,gprime_ext,marker=marker,c=color)
    #ax1.plot(localdens,zz)
    #ax1.axhline(zz[zpyci])
    return gprime_ext

def saltBudget(fname,xval,fig,ax1,title="",color="blue",marker="o"):
    data = timeSeries(fname)
    if "tcdw" not in data:
        return 0
    for k in data.keys():
        if k != "ts":
            try:
                data[k] = np.nanmean(data[k][data["ts"]>2])
            except:
                data[k]=np.nan
    if "gprime" not in data.keys():
        data["gprime"] = np.nan
    variables = grabMatVars(fname,("Xeast","Xwest","Yicefront"))
    yice = float(variables["Yicefront"])
    shelf_width = float(variables["Xeast"])-float(variables["Xwest"])
    shiflx = -data["shiflx"]
    xval,shiflx = FStheory(fname,xval)
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))

    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    minvel = (ds.UVEL*(ds.THETA+1.8)).mean(dim="time",skipna=True).min(dim="XG").min(dim="XC").min(dim="Z")#.argmin(dim="YC")
    ht = ds.UVEL.values*(ds.THETA.values+1.8)
    ht = np.nanmean(ht,axis=0)
    ht = np.nanmin(ht,axis=0)
    ht = np.mean(ht,axis=1)
    minvel=ht

    inputname = "/".join(fname.split("/")[:-1])
    icem = icemask(fname,ds)
    variables = grabMatVars(fname,("Hshelf","Xeast","Xwest","randtopog_height","Yicefront","Zcdw_pt_South","tEast","sEast","icedraft","zz","h","yy","xx"))
    icedraft = np.asarray(variables["icedraft"])
    h = np.asarray(variables["h"])
    meltline = ds.SHIfwFlx.mean(dim="time").mean(dim="XC").values
    meltline = minvel
    weighted = np.sum(meltline*np.asarray(range(len(meltline))))/np.sum(range(len(meltline)))
    #iceline = np.mean(icedraft,axis=0)
    #print("velmin: ",ds.UVEL.mean(dim="time",skipna=True).min(dim="XG").min(dim="Z").argmin(dim="YC").values)
    #plt.plot(iceline)
    #ax2 = plt.gca().twinx()
    #plt.title(fname)
    #ax2.plot(meltline)
    #plt.show()
    #print(minvel)
    #
    tNorth = np.asarray(variables["tEast"])[-1,:]
    sNorth = np.asarray(variables["sEast"])[-1,:]
    rho0 = dens(sNorth,tNorth,0)
    zgl = np.nanmin(icedraft)

    max_height = variables["Zcdw_pt_South"][0][0]
    tcline_height = (max_height-75)/2.0+75
    localdens = dens(sNorth,tNorth,abs(tcline_height))
    d = localdens
    Z = np.asarray(list(ds.Z))
    zz = np.asarray(variables["zz"])[0]
    if np.sum(abs(tNorth-tNorth[0])>0.2)>0:#and np.sum(t>0.5)>0:
        mldi = np.where(abs(tNorth-tNorth[0])>0.2)[0][0]
        #cdwi = np.where(t>0)[0][0]
        rho_1 = np.nanmean(d[:mldi])
        rho_2 = np.nanmean(d[mldi:min(mldi*2,len(d)-1)])
        #gprime_ext = np.nanmean(np.abs(np.diff(d)/np.diff(zz)))
        gprime_ext = 9.8*(rho_2-rho_1)/np.mean((rho_1,rho_2))

    zz = np.asarray(variables["zz"])[0]
    zpyci = np.argmin(np.abs(np.abs(zz)-abs(max_height)))
    gprime_ext = 9.8*(np.mean(localdens[zpyci:min(zpyci*2,len(localdens)-1)])-np.mean(localdens[:zpyci]))/np.mean(localdens[:min(zpyci*2,len(localdens)-1)])
    f = 1.3*10**-4
    shelf_width = float(variables["Xeast"])-float(variables["Xwest"])
    #ax1.scatter(float(data["scdw"]-data["ssurf"])/float(data["tcdw"]-data["tsurf"]),data["froude"],marker=marker,c=color)
    ax1.scatter(shelf_width,float(data["tcdw"]-data["tsurf"])**2,marker=marker,c=color)
    #print(sNorth)
    #ax1.scatter(float(data["scdw"]-data["ssurf"])/float(data["tcdw"]-data["tsurf"]),np.std(sNorth),marker=marker,c=color)
    plt.xlabel(r'Width (m)',fontsize=18)
    plt.ylabel(r'$(T_{cdw} - T_{surf})^2$',fontsize=18)
    return gprime_ext




def steadyStateAverageSimple(fname,xval,fig,ax1,title="",color="blue",marker="o"):
    data = timeSeries(fname)
    for k in data.keys():
        if k != "ts":
            try:
                data[k] = np.nanmean(data[k][data["ts"]>5])
            except:
                data[k]=np.nan
    if "gprime" not in data.keys():
        data["gprime"] = np.nan

    glib = GLIBfromFile(matVarsFile(fname))
    aisf = aisfdepth(matVarsFile(fname))
    shiflx = -data["shiflx"]
    if xval == None:
        xval,shiflx = FStheory(fname,xval)
    ax1.scatter(xval,shiflx,c=color,marker=marker)
    return xval



def steadyStateHT(fname,xval,fig,ax1,title="",color="blue",marker="o"):
    data = timeSeries(fname)
    for k in data.keys():
        if k != "ts":
            try:
                data[k] = np.nanmean(data[k][data["ts"]>2])
            except:
                data[k]=np.nan
    if "gprime" not in data.keys():
        data["gprime"] = np.nan
    variables = grabMatVars(fname,("Xeast","Xwest","Yicefront"))
    yice = float(variables["Yicefront"])
    shelf_width = float(variables["Xeast"])-float(variables["Xwest"])
    shiflx = -data["shiflx"]
    xval,shiflx = FStheory(fname,xval)
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))

    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    inputname = "/".join(fname.split("/")[:-1])
    icem = icemask(fname,ds)

    ax1.scatter(xval,-float(data["incavity"])/np.sum(icem),c=color,marker=marker)
    #ax1.scatter(xval,np.sum(icem),c=color,marker=marker)
    return data["incavity"]

def timeSeriesDashboard(fname,label,fig,axises,times=np.array([]),color="yellow"):
    ((ax1,ax2,ax5),(ax3,ax4,ax6)) = axises 
    data = timeSeries(fname)
    starttime = 0
    if len(data["ts"])==0:
        print(fname, "EMPTY")
        return 
    if np.nanmax(data["ts"][~np.isnan(data["theta"])])<9:
        print(fname,np.nanmax(data["ts"][~np.isnan(data["theta"])]))


    ax1.plot(data["ts"][data["ts"]>starttime],data["theta"][data["ts"]>starttime],label=label)
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Potential Temperature")
    ## salt plot
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

    if np.max(data["ts"][~np.isnan(data["shiflx"])])<7.5:
        print(label)
    
    ## kinetic energy plot
    #totalSHIfwFlx = (ds.SHIfwFlx*ds.XC*ds.YC*ds.Z*ds.hFacC).sum(axis=[1,2,3]).values
    bottomtemps = []
    surfacetemps = []

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
                result = result.group(1)
                if result == "0.25000+02":
                    result = "0.25000e+02"
                return float(result)

def barotropic_streamfunction_graph(fname,description,times=np.array([]),res=1):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    #times=getIterNums(fname)
    #ds = open_mdsdataset(fname,prefix=["UVEL","VVEL"],ignore_unknown_vars=True)#,iters=times)
    ds = open_mdsdataset(fname,prefix=["SALT","THETA","UVEL","VVEL","WVEL"],ignore_unknown_vars=True,extra_variables = extra_variables)
    shortname, fpath = outPath(fname) 
    depth = depthFromdZ(ds)
    UFULL = ds.UVEL.values#
    VFULL = ds.VVEL.values
    #cavity_mask = (ds.SHIfwFlx[0].values !=0)
    #for k in tqdm(range(int(ds.UVEL.values.shape[0]*3/4),ds.UVEL.values.shape[0],res)):
    for k in tqdm(range(10)):
        fig,ax1 = plt.subplots()
        U = UFULL[k]
        V = VFULL[k]
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
        m,s = np.nanmedian(bt[ds.YC.values<60000]),np.nanstd(bt[ds.YC.values<60000])
        #vmin,vmax = m-3*s,m+3*s
        vmin,vmax=0,30000
        frame = ax1.pcolormesh(ds.XC.values,ds.YC.values,bt,vmin=vmin,vmax=vmax)
        ax1.contour(ds.XC.values,ds.YC.values,bt,levels=60,colors="black",vmin=vmin,vmax=vmax)
        ax1.quiver(ds.XC.values,ds.YC.values,np.sum(U,axis=0),np.sum(V,axis=0))
        plt.title(str(k))
        cb = plt.colorbar(frame)
        plt.show()



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

def circulationFigure(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,prefix=["THETA","SALT","momKE","SHIfwFlx","VVEL","UVEL","WVEL"],ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    fig,ax1 = plt.subplots()
    mask = np.logical_and(ds.hFacC.values[0]==0,np.sum(ds.hFacC.values,axis=0)!=0)
    mask[1,:]=0
    zs = ds.Z.values
    xs = ds.XG.values
    ys = ds.YG.values
    ts = np.asarray(times/60.0/60.0/24.0/365.0)
    ts = ds.time.values*grabDeltaT(fname)/60.0/60.0/24.0/365.0
    tsnew = np.full_like(ts,0,dtype=float)
    tsnew[:] = ts
    ts = tsnew
    ts = ts/1000000000

    glib = GLIBfromFile(matVarsFile(fname))

    uvel = np.mean(ds.UVEL.values[ts>5],axis=0)
    vvel = np.mean(ds.VVEL.values[ts>5],axis=0)

    theta= np.mean(ds.THETA.values[ts>5],axis=0)
    hfac = ds.hFacC.values
    interfacez = np.zeros(theta.shape[1:])

    interfaceu = np.zeros(theta.shape[1:])
    interfacev = np.zeros(theta.shape[1:])
    bottomz = np.zeros(theta.shape[1:])
    for i in range(theta.shape[1]):
        for j in range(theta.shape[2]):
            zc = np.where(np.diff(np.sign(theta[:,i,j]-0.5)))[0]
            if len(zc)>0:
                interfacez[i,j] = zs[zc[0]]
                interfaceu[i,j] = uvel[zc[0],i,j]
                interfacev[i,j] = vvel[zc[0],i,j]
            else:
                interfacez[i,j] = np.nan
                interfaceu[i,j] = np.nan
                interfacev[i,j] = np.nan
            if np.sum(hfac[:,i,j]!=0)>0:
                bottomz[i,j] = np.nanmin(zs[hfac[:,i,j]!=0])
            else:
                bottomz[i,j] = np.nan

            
    ax1.set_facecolor("#BBAF98")
    xs,ys=xs/1000,ys/1000
    X,Y= np.meshgrid(xs,ys)
    c= plt.pcolormesh(xs,ys,interfacez[:,::-1],cmap=cmocean.cm.deep)

    caxout = plt.colorbar(c,aspect=40,shrink=0.8,ticks=range(-900,-199,175))
    caxout.ax.tick_params(labelsize=18)
    plt.quiver(X[::5,::5],Y[::5,::5],interfaceu[::5,::5],interfacev[::5,::5],color="white")

    plt.gca().tick_params(labelsize=15)
    plt.contour(xs,ys,bottomz[:,::-1],[glib-25],colors=["red"],linewidths=4)
    plt.ylabel(r'x (km)',fontsize=18)
    plt.xlabel(r'y (km)',fontsize=18)
    plt.tight_layout
    plt.show()



def meltmapmovie(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
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
            frame = ax1.pcolormesh(melt,cmap="jet")
            if k==0:
                cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            #cb.remove()
            frame.remove()


def topMap(fname,description,times=np.array([])):
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))

    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    inputname = "/".join(fname.split("/")[:-1])
    bmask = bottomMask(inputname,ds)
    icem = icemask(inputname,ds)
    shortname, fpath = outPath(fname) 
    fig,ax1 = plt.subplots(1,1)
    theta = ds.THETA.values
    X = ds.XC.values
    Y = ds.YC.values
    thetatop = np.full_like(theta[1,1],0)
    count=0
    for k in tqdm(range(int(theta.shape[0]))):
        if times[k]>5:
            count+=1
            thetak = theta[k]
            thetatopc = thetak *icem
            thetatopc[~icem]=np.nan
            #thetatopc=thetak[0]
            if np.sum(~np.isnan(thetatopc)) !=0:
                thetatop += np.nanmean(thetatopc,axis=0)
    thetatop = thetatop/count
    mask = np.logical_and(ds.hFacC.values[0]==0,np.sum(ds.hFacC.values,axis=0)!=0)
    thetatop[~mask]=np.nan
    
    cax = ax1.pcolormesh(X/1000,Y/1000,thetatop,cmap=cmocean.cm.thermal,vmin=-0.8,vmax=0.2)
    cb = plt.colorbar(cax,ax=ax1,pad=0.05)

    plt.xlabel("x (km)",fontsize=18)
    plt.ylabel("y (km)",fontsize=18)
    plt.xlim(50,350)
    plt.ylim(0,150)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    #ax1.set_title("Width: 100",fontsize=18)
    plt.tight_layout
    plt.show()


    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))

    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    inputname = "/".join(fname.split("/")[:-1])
    bmask = bottomMask(inputname,ds)
    icem = icemask(inputname,ds)
    shortname, fpath = outPath(fname) 
    fig,(ax1,ax2) = plt.subplots(1,2)
    theta = ds.THETA.values
    uvel = ds.UVEL.values
    vvel = ds.VVEL.values
    X,Y = np.meshgrid(range(uvel.shape[3]),range(uvel.shape[2]))
    thetabot,ubot,vbot = np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0)
    thetatop,utop,vtop = np.full_like(theta[1,1],0),np.full_like(theta[1,1],0),np.full_like(theta[1,1],0)
    count=0
    for k in tqdm(range(int(ds.UVEL.values.shape[0]))):
        if times[k]>5:
            count+=1
            thetak = theta[k]
            thetabotc = thetak *bmask
            thetabotc[~bmask] = np.nan
            thetabot += np.nanmean(thetabotc,axis=0)
            ubot += np.nanmean(uvel[k] * bmask,axis=0)
            vbot +=  np.nanmean(vvel[k] * bmask,axis=0)
            thetatopc = thetak *icem
            thetatopc[~bmask] = np.nan
            thetatop += np.nanmean(thetatopc,axis=0)
            utop +=  np.nanmean(uvel[k] * icem,axis=0)
            vtop +=  np.nanmean(vvel[k] * icem,axis=0)
    thetabot = thetabot/count
    thetatop = thetatop/count
    ax1.pcolormesh(X,Y,thetabot,vmin=-2,vmax=2,cmap="jet")

    ax2.pcolormesh(X,Y,thetatop,vmin=-2,vmax=2,cmap="jet")
    
    #ax1.quiver(X,Y,ubot,vbot,scale=0.05)
    #ax2.quiver(X,Y,utop,vtop,scale=0.05)

    ax1.set_title("Bottom")
    ax2.set_title("Top")

    if k==1:
        cb = plt.colorbar(frame,ax=ax2,pad=0.05)



    plt.show()

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
    Znew = np.zeros(ds.THETA.shape[1:])
    for l in range(Znew.shape[0]):
        Znew[l,:,:]=ds.Z.values[l]
    bottom_dist = Znew-(-h.T)
    return np.logical_and(np.logical_and(bottom_dist < 50,bottom_dist>0),ds.hFacC.values>0)

def icemask(fname,ds,thresh=10):
    vals = grabMatVars(fname,("h","icedraft"))
    icedraft = np.abs(np.asarray(vals["icedraft"]))
    Znew = np.zeros(ds.THETA.shape[1:])
    for l in range(Znew.shape[0]):
        Znew[l,:,:]=ds.Z.values[l]
    bottom_dist = Znew-(-icedraft.T)

    return np.logical_and(np.logical_and(bottom_dist > -20,bottom_dist<0),ds.hFacC.values>0)


def mixedLayerQuant(ds,fname):
    THETA = ds.THETA.values
    SALT = ds.SALT.values
    UVEL = ds.UVEL.values
    VVEL = ds.VVEL.values
    vals = grabMatVars(fname,("icedraft","Zcdw_pt_South"))
    icedraft = np.abs(np.asarray(vals["icedraft"]))
    max_height = float(vals["Zcdw_pt_South"][0][0])
    tcline_height = (max_height-75)/2.0+75
    idmask=icedraft>0
    gprimes = []
    ssurf = []
    scdw = []
    tsurf = []
    tcdw = []

    Zfull = np.asarray(list(ds.Z))
    froudes = []
    for t_index in tqdm(range(THETA.shape[0])):
        gprimes_t = []
        ssurf_t = []
        scdw_t = []
        tsurf_t = []
        tcdw_t = []
        froude_map = np.full_like(THETA[0,0,:,:],np.nan)
        for y_index in range(THETA.shape[2]):
            for x_index in range(THETA.shape[3]):
                tcast = THETA[t_index,:,y_index,x_index]
                scast = SALT[t_index,:,y_index,x_index]
                ucast = UVEL[t_index,:,y_index,x_index]
                vcast = VVEL[t_index,:,y_index,x_index]
                if (np.sum(abs(scast)>0.1)>2) and idmask[x_index,y_index]:
                    t = tcast[scast>0.1]
                    Z = Zfull[scast>0.1]
                    u = ucast[scast>0.1]
                    v = vcast[scast>0.1]
                    s = scast[scast>0.1]
                    #t = tcast
                    #s = scast
                    if np.sum(abs(t-t[0])>0.1)>0:#and np.sum(t>0.5)>0:
                        mldi = np.where(abs(t-t[0])>0.1)[0][0]
                        d = dens(s,t,Z[mldi])
                        #cdwi = np.where(t>0)[0][0]
                        rho_1 = np.nanmean(d[:mldi])
                        rho_2 = np.nanmean(d[mldi:max(mldi*2,len(d)-1)])
                        gprimes_t.append(9.8*(rho_2-rho_1)/np.mean((rho_1,rho_2)))
                        ssurf_t.append(np.nanmean(s[:mldi]))
                        scdw_t.append(np.nanmean(s[mldi:max(mldi*2,len(s)-1)]))
                        tsurf_t.append(np.nanmean(t[:mldi]))
                        tcdw_t.append(np.nanmean(t[mldi:max(mldi*2,len(t)-1)]))
                        velmagcdw = np.nanmean(np.sqrt(u[mldi:max(mldi*2,len(t)-1)]**2+v[mldi:max(mldi*2,len(t)-1)]**2))
                        velmagsurf = np.nanmean(np.sqrt(u[:mldi]**2+v[:mldi]**2))
                        
                        froude_map[y_index,x_index]= np.nanmean(np.diff(d)/np.diff(Z))#Z[0]-Z[-1]#np.sqrt((velmagcdw/np.sqrt(gprimes_t[-1]*np.ptp(Z[mldi:max(mldi*2,len(t)-1)])))**2 + (velmagsurf/np.sqrt(gprimes_t[-1]*np.ptp(Z[:mldi+1])))**2)
                        #froude_map[y_index,x_index]= np.ptp(Z)

        #print(fname)
        #plt.imshow(froude_map,vmin=0,vmax=1,cmap="jet")
        #plt.pcolormesh(ds.XC.values,ds.YC.values,froude_map,cmap="jet")
        #plt.colorbar()
        #plt.show()
        #froudes.append(np.nanmean(froude_map))

        gprimes.append(np.nanmean(gprimes_t))
        ssurf.append(np.nanmean(ssurf_t))
        scdw.append(np.nanmean(scdw_t))
        tsurf.append(np.nanmean(tsurf_t))
        tcdw.append(np.nanmean(tcdw_t))
    #plt.imshow(gprimemap)
    #plt.show()
    return gprimes,ssurf,scdw,tsurf,tcdw,froudes#,froudesurf,froudecdw


    
def outPath(resultspath):
    nameparts = resultspath.split("/")
    shortname = nameparts[-3] + "|" + nameparts[-2]
    fpath = "/".join((nameparts[:-5]+["pics"] + [shortname]))
    return shortname, fpath
    
def meltMapAverage(fname,description,quant="THETA",res=1,dim="zonal",ax1=None):
    if not ax1:
        fig,ax1 = plt.subplots(figsize=(10,8))
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    ds = open_mdsdataset(fname,ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)
    mask = np.logical_and(ds.hFacC.values[0]==0,np.sum(ds.hFacC.values,axis=0)!=0)
    times = times*grabDeltaT(fname)/60.0/60.0/24.0/365.0
    melt= -np.nanmean(ds.SHIfwFlx.values[times>2.5],axis=0)*(60*60*24*365)*(1/920.0)
    xs = ds.XG.values/1000
    ys = ds.YG.values/1000
    newcmap = cmocean.tools.crop(cmocean.cm.balance, 0, 45, 0)

    vals = grabMatVars(fname,("h","icedraft"))
    icedraft = np.abs(np.asarray(vals["icedraft"]))
    h = np.abs(np.asarray(vals["h"]))
    icedraft = np.logical_and(icedraft!=0,icedraft!=h)
    melt[~icedraft.T]=np.nan
    
    im = ax1.pcolormesh(xs,ys,melt,cmap=newcmap,vmin=0,vmax=25)

    bound = np.argwhere(~np.isnan(melt))

    #ax2.set_xlim((140,260))
    #ax2.set_ylim((0,160))
    pad = 10
    ax1.set_xlim(xs[min(bound[:, 1])]-pad, xs[max(bound[:, 1])]+pad)
    ax1.set_ylim(ys[min(bound[:, 0])]-pad, ys[max(bound[:, 0])]+pad)

    ax1.set_xlabel("x (km)",fontsize=18)
    ax1.set_ylabel("y (km)",fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    caxout = inset_axes(
        ax1,
        width="2%",  # width: 5% of parent_bbox width
        height="40%",  # height: 50%
        loc="lower left",
        bbox_to_anchor=(1.0, 0.30, 1, 1),
        bbox_transform=ax1.transAxes,
        borderpad=1,
    )
    caxout.tick_params(labelsize=18)
    plt.colorbar(im,cax=caxout,ticks=[0,20])
    caxout.set_ylabel('m/yr', rotation=0,fontsize=18)

    ax1.set_ylabel('Y (km)',fontsize=18)
    ax1.set_xlabel('X (km)',fontsize=18)

    #caxout = plt.colorbar(im,ax=ax1, aspect=40,shrink=0.4,ticks=range(0,41,10))
    #caxout.ax.tick_params(labelsize=18)
    #plt.show()


def crossSectionAverage(fname,description,quant="THETA",res=1,dim="zonal",ax1=None,show=True):
    if not ax1:
        fig,ax1 = plt.subplots(figsize=(10,8))
    extra_variables = dict( SHIfwFlx = dict(dims=["k","j","i"], attrs=dict(standard_name="Shelf Fresh Water Flux", units="kg/m^3")))
    times=getIterNums(fname)
    if quant!="DENS":
        ds = open_mdsdataset(fname,prefix=["THETA"],ignore_unknown_vars=True,extra_variables = extra_variables,iters=times)

        if dim == "zonal":
            #zonal_average = ds.where(ds.hFacC==1).mean(dim="XC",skipna=True)
            zonal_average = ds.sel(XC=200*10**3,method="nearest",drop=True)
            ys = zonal_average.YC.values
        if dim == "meridional":
            zonal_average = ds.isel(YC=50)
            ys = zonal_average.XC.values
        zvals = zonal_average[quant].values
        #zvals[:,zonal_average.hFacC.values!=1]=np.nan
        #zvals[zvals==0]=np.nan
        m = np.nanmedian(zvals)
        s = np.nanstd(zvals)
        tmin, tmax = m-2*s,m+s*2
        shortname, fpath = outPath(fname) 
        #plt.hist(zvals[:,:,:].flatten())
        #plt.show()
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

    times = times*grabDeltaT(fname)/60.0/60.0/24.0/365.0
    zvals=np.nanmean(zvals[times>2.5],axis=0)
    newys = []
    for i in range(len(ys)-1):
        newys.append(np.linspace(ys[i],ys[i+1],10))
    newzs = []
    for i in range(len(zs)-1):
        newzs.append(np.linspace(zs[i],zs[i+1],10))
    newys = np.concatenate(newys)
    newzs = np.concatenate(newzs)
    newzvals = np.empty((np.asarray(zvals.shape)-1)*10)
    facs = zonal_average.hFacC.values
    for i in range(zvals.shape[0]):
        for j in range(zvals.shape[1]):
            if facs[i,j]==1 or (facs[i,j]!=1 and facs[max(i-1,0),j]==1):
                newzvals[i*10:i*10+int(facs[i,j]*10),j*10:j*10+10]=zvals[i,j]
                newzvals[i*10+int(facs[i,j]*10):i*10+10,j*10:j*10+10]=np.nan
            elif facs[i,j]==0:
                newzvals[i*10:i*10+10,j*10:j*10+10]=np.nan
            elif facs[i,j]!=1 and facs[min(i+1,zvals.shape[0]-1),j]==1:
                newzvals[i*10:i*10+10-int(facs[i,j]*10),j*10:j*10+10]=np.nan
                newzvals[i*10+10-int(facs[i,j]*10):i*10+10,j*10:j*10+10]=zvals[i,j]
            #elif facs[i,j]==0:
                #newzvals[i,j]=np.nan
    #c=ax1.pcolormesh(newys/1000,newzs/1000,newzvals,cmap=cmocean.cm.thermal)
    #for j in range(zvals.shape[1]):
        #nonzero = np.nonzero(facs[:,j])
        #if len(nonzero[0])>0:
            #i = int(nonzero[0][0])
        ##else:
            #print("zoop")
            #i = np.where(-newzs>1500)[0][0]
        #newzvals[i*10:,j*10:j*10+10]=np.nan
        #newzvals[:i*10,j*10:j*10+10]=1
    c=ax1.pcolormesh(newys/1000,newzs/1000,newzvals,cmap=cmocean.cm.thermal,zorder=5,vmin=-2,vmax=1)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    variables = grabMatVars(fname,("Xeast","Xwest","Yicefront"))
    yice = float(variables["Yicefront"])/1000
    ices = volume(fname)
    C=(-.205)-yice*ices
    y = np.array([[0.05,0], [0.05,0.05*ices+C], [yice + 0.8,-.205-0.1], [yice + 0.8,0], [0.05,0]])
    p = Polygon(y, facecolor = '#8addf9',zorder=4)
    ax1.add_patch(p)

        
    ax1.set_facecolor("#BBAF98")
    ax1.set_ylim(-1.5,0.0051)
    ax1.set_xlim(-0.5,299)
    caxout = inset_axes(
        ax1,
        width="2%",  # width: 5% of parent_bbox width
        height="40%",  # height: 50%
        loc="lower left",
        bbox_to_anchor=(1.0, 0.30, 1, 1),
        bbox_transform=ax1.transAxes,
        borderpad=1,
    )
    caxout.tick_params(labelsize=18)
    plt.colorbar(c,cax=caxout,ticks=[-2,0,1])
    caxout.set_ylabel('$^\circ$C', rotation=0,fontsize=18)

    ax1.set_ylabel('Depth (km)',fontsize=18)
    ax1.set_xlabel('Cross Shelf Distance (km)',fontsize=18)



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
        if dim == "zonal":
            #zonal_average = ds.where(ds.hFacC==1).mean(dim="XC",skipna=True)
            zonal_average = ds.isel(XC=75)
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
    with moviewriter.saving(fig, fpath+quant+"|"+dim+".mp4" , dpi=250):
        print(fpath+quant+"|"+dim+".mp4")
        for k in tqdm(range(0,length,res)):
            if quant == "DENS":
                frame = ax1.pcolormesh(ys,zs,gsw.sigma0(zvals[0][k,:,:],zvals[1][k,:,:]),cmap="jet",vmin=27.4,vmax=27.7)
            elif quant == "THETA":
                frame = ax1.pcolormesh(ys,zs,zvals[k,:,:],cmap="jet",vmin=-0.5,vmax=1)
            elif quant == "SALT":
                frame = ax1.pcolormesh(ys,zs,zvals[k,:,:],cmap="jet",vmin=34,vmax=34.8)
            else:
                frame = ax1.imshow(zvals[k,:,:],cmap="jet")
            cb = plt.colorbar(frame)
            moviewriter.grab_frame()
            cb.remove()
            frame.remove()
    plt.close()

def bottomAnim(fname,description,times=np.array([]),quant="THETA",res=5):
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
            frame = ax1.pcolormesh(ds.XC.values,ds.YC.values,znew,cmap="jet",vmin=-2,vmax=1)
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

def legendFunction(runsdict):
    conversion={"d0":"cavityd0","d-200":"cavityd-200","slope200":"cavityd200","slope375":"cavityd375"}
    for k in runsdict.keys():
        for l in range(len(runsdict[k]["specialstring"])):
            if not runsdict[k]["specialstring"][0]:
                ss = k.split("-")[0]
            else:
                ss = runsdict[k]["specialstring"][l]
            if ss in conversion.keys():
                ss = conversion[ss]
            plt.scatter(1,1,marker=runsdict[k]["marker"][l],color=runsdict[k]["color"][l],label=ss)
    plt.legend()
    plt.show()
            

def folderMapTimeSeries(runsdict,save=True):
    fig,axises = plt.subplots(2,3,figsize=(8,7))
    for k in runsdict.keys():
        for f in glob.glob(str("/home/garrett/Projects/MITgcm_ISC/experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                if key and key in f:
                    try:
                        timeSeriesDashboard(f+"/results",key+f[-6:-10],fig,axises,color=runsdict[k]["color"][l])
                    except:
                        print("yeesh")
                else:
                    try:
                        timeSeriesDashboard(f+"/results","",fig,axises,color=runsdict[k]["color"][l])
                    except:
                        print("yeesh")
    axises[0][0].legend()

    plt.show()

def fastExplosionCheck(runsdict,save=True):
    fig,axises = plt.subplots(2,3,figsize=(8,7))
    count = [0,0]
    for k in tqdm(runsdict.keys()):
        for f in glob.glob(str("/home/garrett/Projects/MITgcm_ISC/experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                if key and key in f:
                    
                    try:
                        with open(f+"/results/STDOUT.0000") as myfile:
                            if 'NaN' in myfile.read():
                                print(f, "oh no a NaN!")
                            else:
                                count[0]=count[0]+1
                            count[1]=count[1]+1
                    except:
                        print("no STDOUT")
    print(count)


def folderMap(runsdict,save=True):
    fig,axises = plt.subplots(1,1,figsize=(8,7))
    xs,ys,eyeds = [],[],{}
    for k in runsdict.keys():
        for f in glob.glob(str("/home/garrett/Projects/MITgcm_ISC/experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                if key and key in f:
                    try:
                        x,y=FStheory(f+"/results",None)
                        if ~np.isnan(y):
                            xs.append(x)
                            ys.append(y)
                            eyeds[f+str(l)]=len(xs)-1
                    except:
                        print("yeesh")
                elif not key:
                    try:
                        x,y=FStheory(f+"/results",None)
                        if ~np.isnan(y):
                            xs.append(x)
                            ys.append(y)
                            eyeds[f+str(l)]=len(xs)-1
                    except:
                        print("yeesh")
    xs = np.asarray(([xs])).reshape((-1, 1))
    model = LinearRegression().fit(xs, ys)
    rho0 = 1025
    rhoi = 910
    Cp = 4186
    If = 334000
    C = model.coef_
    #W0 = (rho0*Cp)/(rhoi*If*C)
    W0 =  100000#(rho0*Cp)/(rhoi*If*C)
    alpha =  C/((rho0*Cp)/(rhoi*If*W0))
    print("alpha:",alpha)
    oldxs=xs
    #xs= xs*(rho0*Cp)/(rhoi*If*325000)
    xs=model.predict(xs)
    plt.text(.05, .95, '$r^2=$'+str(round(model.score(oldxs,ys),2)), ha='left', va='top', transform=plt.gca().transAxes,fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.xlabel("Melt rate predicted from theory (m/yr)",fontsize=18)
    plt.ylabel("Diagnosed melt rate from simulation (m/yr)",fontsize=18)
    plt.plot(range(0,25),range(0,25),linestyle="dashed")
    plt.xlim(0,20)
    plt.ylim(0,20)
    for k in runsdict.keys():
        for f in glob.glob(str("/home/garrett/Projects/MITgcm_ISC/experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                if key and key in f:
                    if f+str(l) in eyeds.keys():
                        steadyStateAverageSimple(f+"/results",xs[eyeds[f+str(l)]]*(60*60*24*365),fig,axises,color=runsdict[k]["color"][l],marker=runsdict[k]["marker"][l],title=runsdict[k]["description"][0])
                elif not key:
                    try:
                        if f+str(l) in eyeds.keys():
                            steadyStateAverageSimple(f+"/results",xs[eyeds[f+str(l)]]*(60*60*24*365),fig,axises,color=runsdict[k]["color"][l],marker=runsdict[k]["marker"][l],title=runsdict[k]["description"][0])
                    except:
                        print("yeesh")
            if save:
                plt.savefig("/home/garrett/Projects/HUB/paperfigures/"+k+".png")

def folderMapGeneric(func,runsdict,save=True):
    fig,axises = plt.subplots(1,1,figsize=(8,7))
    for k in runsdict.keys():
        for f in glob.glob(str("/home/garrett/Projects/MITgcm_ISC/experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                if key and key in f:
                    try:
                        func(f+"/results",None,fig,axises,color=runsdict[k]["color"][l],marker=runsdict[k]["marker"][l],title=runsdict[k]["description"][0])
                    except:
                        print("yeesh")
                elif not key:
                    #try:
                    func(f+"/results",None,fig,axises,color=runsdict[k]["color"][l],marker=runsdict[k]["marker"][l],title=runsdict[k]["description"][0])
                    #except:
                        #print("yeesh")

            if save:
                plt.savefig("/home/garrett/Projects/HUB/paperfigures/"+k+".png")

def folderMapMoreGeneric(func,runsdict):
    for k in runsdict.keys():
        for f in glob.glob(str("/home/garrett/Projects/MITgcm_ISC/experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                print(f)
                key=runsdict[k]["specialstring"][l]
                if "at125" in f:
                    if key and key in f:
                        #try:
                        func(f+"/results",key)
                        #except:
                        #print("yeesh")
                    elif not key:
                        #try:
                        func(f+"/results",key)
                        #except:
                            #print("yeesh")


def folderMapRefresh(runsdict,save=True):
    fig,axises = plt.subplots(1,1,figsize=(8,7))
    for k in runsdict.keys():
        for f in glob.glob(str("/home/garrett/Projects/MITgcm_ISC/experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                if key and key in f :
                    try:
                        timeSeries(f+"/results",True)
                    except:
                        print("yeesh")
                elif not key:
                    try:
                        if f+str(l) in eyeds.keys():
                            timeSeries(f+"/results",True)
                    except:
                        print("yeesh")

            if save:
                plt.savefig("/home/garrett/Projects/HUB/paperfigures/"+k+".png")
       
runsdict = {\
        "icefront-64":{"specialstring":['y100','y250'], "marker":["^","^"] ,"color":["orange","purple"],"description":["Different ice front distances"]},\
        "shelfdepth-16":{"specialstring":['d500','d600','d700'], "marker":["x","x","x","x"] ,"color":["purple","gray","green","red"],"description":["Different shelf depths"]},\
        "widthexp-GLIB-explore-32":{"specialstring":['w50','w100','w250'], "marker":["o","o","o"] ,"color":["orange","purple","gray"],"description":["Different shelf widths" ]},\
        "slope-22":{"specialstring":["s300","s150","s0"], "marker":["D","D","D"] ,"color":["red","green","blue"],"description":["more steep slope"]},\
        "reference":{"specialstring":["ref"], "marker":["p"] ,"color":["black"],"description":["more steep slope"]}\
        }

#generateRunsTable(runsdict)
#legendFunction(runsdict)
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-32/at0w250/results/","")

#fastExplosionCheck(runsdict)
#folderMapRefresh(runsdict)
#folderMapTimeSeries(runsdict)
folderMap(runsdict)
#folderMapGeneric(gprimeWidth,runsdict)
plt.show()
#plt.show()
#crossSectionAverage("/home/garrett/Projects/MITgcm_ISC/experiments/reference/at125/results","Reference")
#plt.show()
#circulationFigure("/home/garrett/Projects/MITgcm_ISC/experiments/reference/at125/results","Reference")
#plt.show()

######FIX ITERNUMS PROBLEM
#folderMapGeneric(steadyStateAverageSimple,runsdict)
#folderMapGeneric(gprimeWidth,runsdict)
#folderMapGeneric(saltBudget,runsdict)
#plt.show()
#folderMapGeneric(steadyStateAverageSimple,runsdict)
#plt.show()

