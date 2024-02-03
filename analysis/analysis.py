from datainput import *
from matlabglib import GLIBfromFile
import numpy as np
from jmd95 import dens
from scipy.integrate import quad

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

def slope(fname,method="grad"):
    if method == "param":
        variables = grabMatVars(fname,('Zcdw_pt_shelf','icedraft','tEast','zz','yy',"xx","Yicefront"))
        y = np.asarray(variables["Yicefront"])[0][0]
        icedraft = np.asarray(variables["icedraft"])
        zgl = np.nanmin(icedraft)-np.nanmax(icedraft[icedraft!=0])
        return (abs(zgl)-200)/y
    if method == "grad":
        variables = grabMatVars(fname,('icedraft',"h","YY","xx","yy","Yicefront"))
        icedraft = np.asarray(variables["icedraft"])
        h = np.asarray(variables["h"])
        YY = np.asarray(variables["YY"])
        yy = np.asarray(variables["yy"])[0]
        xx = np.asarray(variables["xx"])[0]
        diff = np.abs(h)-np.abs(icedraft)
        grad = np.gradient(icedraft)
        grad[0] = (grad[0]/np.gradient(xx)[10])**2
        grad[1] = (grad[1]/np.gradient(yy)[10])**2
        #grad = np.sqrt(grad[0] + grad[1])
        grad = np.sqrt(np.sum(grad,axis=0))
        return np.nanmedian(grad[np.logical_and(icedraft!=0,diff!=0)])#np.mean(diff[np.logical_and(icedraft!=0,diff!=0)]) #+ abs(zglib-zgl)/y

def FStheory(fname,xval):
    data = timeSeries(fname)
    glib = GLIBfromFile(matVarsFile(fname))
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
    glibxval = intTemp(glib,fname)
    ices = slope(fname)
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
    rho_1i = np.logical_and(zz<zz[zpyci],zz>zz[zpyci]-30)
    rho_2i = np.logical_and(zz<zz[zpyci]+30,zz>zz[zpyci])

    gprime_ext = 9.8*(np.mean(localdens[rho_1i])-np.mean(localdens[rho_2i]))/np.mean(localdens[np.logical_or(rho_1i,rho_2i)])

    deltaH = -(abs(tcline_height)- abs(glib))
    if "reference" in fname and "at125" in fname:
        print(tcline_height)

    f = 1.3*10**-4
    rho0 = 1025
    rhoi = 910
    Cp = 4186
    If = 334000
    return (glibxval)*deltaH*(gprime_ext)/(f)*ices,-data["shiflx"]/(60*60*24*365)

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


