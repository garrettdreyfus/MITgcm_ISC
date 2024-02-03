%%%
%%% calcOverturning_pt_Aocean.m
%%%
%%% Calculates the overturning circulation, calculated using the MITgcm 
%%% 'layers' package. Writes the result to a .mat file.
%%%
%%% expdir - Base directory containing experiment
%%% expname - Name of experiment
%%% prodir - Directory in which to write output .mat file
%%%


function calcOverturning_pt_Aocean (expdir,expname,prodir)

  %%% Load experiment
  loadexp;
  disp(layers_bounds)
  %%% Density bins for MOC calculation  
  ptlevs = layers_bounds(:,2);
  Npt = length(ptlevs)-1;

  %%% Frequency of diagnostic output - should match that specified in
  %%% data.diagnostics.
  dumpFreq = abs(diag_frequency(1));
  nDumps = round(nTimeSteps*deltaT/dumpFreq);
  dumpIters = round((1:nDumps)*dumpFreq/deltaT);
  dumpIters = dumpIters(dumpIters > nIter0);

  %%% Create a finer vertical grid
  ffac = 5;
  Nrf = ffac*Nr;
  delRf = zeros(1,Nrf); 
  for n=1:Nr
    for m=1:ffac
      delRf((n-1)*ffac+m) = delR(n)/ffac;
    end
  end
  zz = - cumsum((delR + [0 delR(1:Nr-1)])/2);
  zz_f = - cumsum((delRf + [0 delRf(1:Nrf-1)])/2);

  %%% Partial cell heights on fine grid
  hFacS_f = zeros(Nx,Ny,Nrf);
  for k=1:Nr
    hFacS_f(:,:,ffac*(k-1)+1:ffac*k) = hFacS(:,:,k*ones(1,ffac));              
  end

  %%% Grid of actual vertical positions, accounting for partial cells
  ZZ = zeros(Nx,Ny,Nr);
  ZZ_f = zeros(Nx,Ny,Nrf);
  DZ = zeros(Nx,Ny,Nr);
  DZ_f = zeros(Nx,Ny,Nrf);
  PP = zeros(Nx,Ny,Nr);
  ZZ(:,:,1) = - delR(1)*hFacS(:,:,1)/2;
  for k=2:Nr
    ZZ(:,:,k) = ZZ(:,:,k-1) - 0.5*delR(k-1)*hFacS(:,:,k-1) - 0.5*delR(k)*hFacS(:,:,k);
  end       
  ZZ_f(:,:,1) = - delRf(1)*hFacS_f(:,:,1)/2;
  for k=2:Nrf 
    ZZ_f(:,:,k) = ZZ_f(:,:,k-1) - 0.5*delRf(k-1)*hFacS_f(:,:,k-1) - 0.5*delRf(k)*hFacS_f(:,:,k);      
  end
  for k=1:Nr
    DZ(:,:,k) = delR(k);
  end   
  for k=1:Nrf
    DZ_f(:,:,k) = delRf(k);
  end   
  for k=1:Nr
    PP(:,:,k) = -delR(k);
  end   

  %%% Matrices for vertical interpolation  
  k_p = zeros(Nx,Ny,Nrf);
  k_n = zeros(Nx,Ny,Nrf);
  w_n = zeros(Nx,Ny,Nrf);
  w_p = zeros(Nx,Ny,Nrf);
  is_wet_col = zeros(Nx,Ny);
  for i=1:Nx
    for j=1:Ny

      %%% Indices of the lowest cells
      kmax = sum(squeeze(hFacS(i,j,:))~=0);
      kmax_f = ffac*kmax;
      is_wet_col(i,j) = (kmax~=0);

      for k=1:Nrf

        %%% Previous and next interpolation indices
        k_p(i,j,k) = ceil(k/ffac-0.5);
        k_n(i,j,k) = k_p(i,j,k) + 1;

        %%% Fine grid cell is above highest coarse grid cell, so fine grid
        %%% gamma will just be set equal to uppermost coarse grid gamma
        if (k_p(i,j,k) <= 0)

          k_p(i,j,k) = 1;
          w_p(i,j,k) = 0;
          w_n(i,j,k) = 1;

        else

          %%% Fine grid cell is below lowest coarse grid cell, so fine grid
          %%% gamma will just be set equal to lowermost coarse grid gamma
          if (k_n(i,j,k) > kmax)

            k_n(i,j,k) = kmax;
            w_n(i,j,k) = 0;
            w_p(i,j,k) = 1;

          %%% Otherwise set weights to interpolate linearly between neighboring
          %%% coarse-grid gammas
          else

            w_p(i,j,k) = (ZZ(i,j,k_n(i,j,k))-ZZ_f(i,j,k))./(ZZ(i,j,k_n(i,j,k))-ZZ(i,j,k_p(i,j,k)));
            w_n(i,j,k) = 1 - w_p(i,j,k);

          end

        end

      end

    end
  end

  %%% Calculate time-averaged isopycnal flux, density and velocity 
  load([prodir expname '_tavg_5yrs.mat'],'LaVH2TH','LaHs2TH','VVEL','THETA');
  vflux_tavg = LaVH2TH;  
  h_pt_tavg = LaHs2TH;   
  vvel_tavg = VVEL;  
  pt_tavg = THETA;  

  %%% Interpolate potential temperature to v-gridpoints  
  pt_v = NaN*pt_tavg;
  pt_v(:,2:Ny,:) = 0.5* (pt_tavg(:,1:Ny-1,:) + pt_tavg(:,2:Ny,:));    

  %%% Interpolate onto a finer grid         
  vvel_f = zeros(Nx,Ny,Nrf);
  pt_f = NaN*zeros(Nx,Ny,Nrf);
  if (ffac == 1)

    %%% Shortcut if fine grid resolution = coarse grid resolution
    vvel_f = vvel_tavg;        
    pt_f = pt_v;

  else   

    %%% Velocity uniform throughout each coarse grid cell to preserve
    %%% mass conservation
    for k=1:Nr
      vvel_f(:,:,ffac*(k-1)+1:ffac*k) = vvel_tavg(:,:,k*ones(1,ffac));          
    end

    %%% Linearly interpolate density
    for i=1:Nx
      for j=1:Ny %%% Restrict to wet grid cells  
        if (is_wet_col(i,j))
          pt_f(i,j,:) = w_p(i,j,:).*pt_v(i,j,squeeze(k_p(i,j,:))) + w_n(i,j,:).*pt_v(i,j,squeeze(k_n(i,j,:)));
        end
      end
    end

  end            

  %%% Calculate mean fluxes within mean density surfaces
  vflux_m = 0*vflux_tavg;
  vdz = vvel_f.*hFacS_f.*DZ_f;
  vflux_m(:,:,Npt) = vflux_m(:,:,Npt) + sum(vdz.*(pt_f>ptlevs(Npt)),3);
  vflux_m(:,:,1) = vflux_m(:,:,1) + sum(vdz.*(pt_f<=ptlevs(2)),3);
  for m=2:Npt-1
    vflux_m(:,:,m) = vflux_m(:,:,m) + sum(vdz.*((pt_f>ptlevs(m)) & (pt_f<=ptlevs(m+1))),3);
  end   

  %%% Zonally integrate meridional fluxes
  vflux_xint = zeros(Ny,Npt);
  vflux_m_xint = zeros(Ny,Npt);
  for i=1:Nx
    vflux_xint = vflux_xint + delX(i)*squeeze(vflux_tavg(i,:,:));
    vflux_m_xint = vflux_m_xint + delX(i)*squeeze(vflux_m(i,:,:));
  end

  %%% Sum fluxes to obtain streamfunction
  psi_pt = zeros(Ny,Npt+1);
  psim_pt = zeros(Ny,Npt+1);
  for m=1:Npt  
    psi_pt(:,m) = sum(vflux_xint(:,m:Npt),2);     
    psim_pt(:,m) = sum(vflux_m_xint(:,m:Npt),2);     
  end
  psi_pt = psi_pt/1e6;
  psim_pt = psim_pt/1e6;
  psie_pt = psi_pt - psim_pt;

  %%% Calculate mean density surface heights
  h_pt_xtavg = squeeze(nanmean(h_pt_tavg));
  z_pt = 0*h_pt_xtavg;
  for m=1:Npt
    z_pt(:,m) = - sum(h_pt_xtavg(:,1:m-1),2); % equivalent to zisop_mean
  end
  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% CALCULATE ISOPYCNAL DEPTHS %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%
    Aocean = zeros(Ny,Nr+1);
    Aisop = zeros(Ny,Npt+1);

    %%% Grid spacing matrices
    DX_xyz = repmat(reshape(delX,[Nx 1 1]),[1 Ny Nr]);
    DZ_xyz = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);
    DX_xypt = repmat(reshape(delX,[Nx 1 1]),[1 Ny Npt]);


    Aocean_xint = squeeze(sum(DX_xyz.*DZ_xyz.*hFacS,1)); %%% Integrate the ocean area in the x-direction, on v-grid
    Aocean(:,1:Nr) = cumsum(Aocean_xint,2,'reverse'); %%% Integrate the ocean area in the z-direction from bottom to top

    Aisop_xint = squeeze(nansum(DX_xypt.*h_pt_tavg,1));
    Aisop(:,2:Npt+1) = cumsum(Aisop_xint,2);

    zzf = -[0 cumsum(delR)];

    figure(30)
    subplot(1,2,1)
    pcolor(yy/1000,zzf,Aocean'/Lx);
    shading interp;colormap('default');colorbar;
    xlabel('y (km)');ylabel('z (m)');
    title('Aocean/Lx (m)')
    subplot(1,2,2)
    pcolor(yy/1000,ptlevs(1:Npt+1),Aisop'/Lx);
    shading interp;colormap('default');colorbar;
    xlabel('y (km)');ylabel('Potential temperature (degC)');
    title('Aisop/Lx (m)')

    %%% The maximum of Aisop should be approximately the same as Aocean
    diff_Aocean_Aisop = (Aocean(:,1)-Aisop(:,end))/Lx; 

    psi_z = zeros(Ny,Nr+1);
    psim_z = zeros(Ny,Nr+1);
    psie_z = zeros(Ny,Nr+1);
    Zisop = zeros(Ny,Npt+1);
    
    for j = 1:Ny
        %%% No ocean here so can't interpolate
        if (Aocean(j,1)==0) 
            continue;
        end
        %%% Truncate Aisop and psi_pt vectors to remove repeated entries at the ends
        %%% of the Aisop vector
        [Aisop_trunc,idx] = unique(Aisop(j,:));
        psi_pt_trunc = psi_pt(j,idx);
        psim_pt_trunc = psim_pt(j,idx);
        %%% Interpolate vertically. Note that Aisop(k) is the total area over
        %%% density bins 1 through k, i.e. it corresponds to the density level
        %%% k+1/2.
        psi_z(j,:) = interp1(Aisop_trunc,psi_pt_trunc,Aocean(j,:),'linear','extrap');
        psim_z(j,:) = interp1(Aisop_trunc,psim_pt_trunc,Aocean(j,:),'linear','extrap');
        
        
        %%% Truncate Aocean and zzf vectors to remove repeated entries at the ends
        %%% of the Aocean vector
        [Aocean_trunc,idx_aocean] = unique(Aocean(j,:));
        zzf_trunc = zzf(idx_aocean);       
        Zisop(j,:) = interp1(Aocean_trunc,zzf_trunc,Aisop(j,:),'linear','extrap');
    end
    
    psie_z = psi_z - psim_z;

    [DD,LL] = meshgrid(ptlevs,yy);

    figure(31)
    subplot(2,2,1)
    pcolor(yy/1000,ptlevs,psi_pt');
    shading interp;colormap('redblue');colorbar;caxis([-5 5]);
    xlabel('y (km)');ylabel('Potential temperature (degC)');
    title('\psi (Sv) in PT space')
    subplot(2,2,2)
    pcolor(yy/1000,zzf,psi_z');
    shading interp;colormap('redblue');colorbar;caxis([-5 5]);
    xlabel('y (km)');ylabel('z (m)');
    title('\psi (Sv), interpolating the streamfunction')
    subplot(2,2,3)
    pcolor(LL/1000,Zisop,psi_pt);
    shading interp;colormap('redblue');colorbar;caxis([-5 5]);
    xlabel('y (km)');ylabel('z (m)');
    title('\psi (Sv), interpolating Zisop')
    subplot(2,2,4)
    pcolor(LL/1000,Zisop,psi_pt);
    shading interp;colormap('redblue');colorbar;caxis(["auto"]);
    hold on;
    contour(LL/1000,Zisop,DD,[-1:1:12],'EdgeColor','w');
    hold off;
    xlabel('y (km)');ylabel('z (m)');
    title('\psi (Sv), interpolating Zisop (with isotherms)')


  %%% Store computed data for later
  save(fullfile(prodir,[expname,'_MOC_pt_Aocean.mat']), ...
    'xx','yy','zz','ptlevs','zzf', ... 
    'pt_tavg', ...  
    'vflux_m','vflux_xint','vflux_m_xint', ...
    'psi_pt','psim_pt','psie_pt', ...
    'psi_z','psim_z','psie_z','DD','LL','Aisop','Aocean','Zisop');
  
end
