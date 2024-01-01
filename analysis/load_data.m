    
    

    nIter = 2815714;
clear oceTAUX oceTAUY
%     if(expname(end-4:end)=='_prod')
    if(false)
        load([prodir expname '_tavg_5yrs.mat'],'THETA','SALT','UVEL','VVEL','VVELTH','UVELTH','ETAN',...
                'SHIfwFlx','oceTAUX','oceTAUY','SHI_TauY','WVEL','WVELTH');
        tt = THETA;
        ss = SALT;
        uu = UVEL;
        vv = VVEL;
        ww = WVEL;
        vt = VVELTH;
        ut = UVELTH;
        wt = WVELTH;
        eta = ETAN;
        %    rho_insitu = RHOAnoma+rhoConst;
        if(useSEAICE)
            load([prodir expname '_tavg_5yrs.mat'],'SIuice','SIvice');
            ui = SIuice;
            vi = SIvice;
        end
    else
        tt = rdmds([exppath,'/results/THETA'],nIter);
        ss = rdmds([exppath,'/results/SALT'],nIter);
        uu = rdmds([exppath,'/results/UVEL'],nIter);
        vv = rdmds([exppath,'/results/VVEL'],nIter);
        ww = rdmds([exppath,'/results/WVEL'],nIter);
        eta = rdmds([exppath,'/results/ETAN'],nIter);
        SHIfwFlx = rdmds([exppath,'/results/SHIfwFlx'],nIter);
        SHI_TauY = rdmds([exppath,'/results/SHI_TauY'],nIter);
        if(useSEAICE)
            ui = rdmds([exppath,'/results/SIuice'],nIter(n));
            vi = rdmds([exppath,'/results/SIvice'],nIter(n));
        end
    end
    
