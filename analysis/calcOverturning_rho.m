%%%
%%% calcOverturning_rho_new.m
%%%
%%% Calculates the overturning circulation, calculated using the MITgcm 
%%% 'layers' package. Writes the result to a .mat file.
%%%
%%% expdir - Base directory containing experiment
%%% expname - Name of experiment
%%% outdir - Directory in which to write output .mat file
%%% tmin - Start time for analysis period, in years
%%% tmax - End time for analysis period, in years
%%%
function calcOverturning_rho (expdir,expname,prodir)
 
%%
  %%% Load experiment
  loadexp;
 
  %%% Density bins for MOC calculation  
%   ptlevs = flip(layers_bounds);
  ptlevs = flip(layers_bounds(:,1));
%   ptlevs = layers_bounds;
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
%   vflux_tavg = readIters(exppath,'LaVH1RHO',dumpIters,deltaT,tmin,tmax,Nx,Ny,Npt);    
%   h_pt_tavg = readIters(exppath,'LaHs1RHO',dumpIters,deltaT,tmin,tmax,Nx,Ny,Npt);    
%   vvel_tavg = readIters(exppath,'VVEL',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);
%   temp_tavg = readIters(exppath,'THETA',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);    
%   salt_tavg = readIters(exppath,'SALT',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);  
%   pressure_tavg = readIters(exppath,'PHIHYD',dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr);  
load([prodir expname '_tavg_5yrs.mat'],'LaVH1RHO','LaHs1RHO','VVEL','THETA','SALT','PHIHYD');
  vflux_tavg = flip(LaVH1RHO,3);  
  h_pt_tavg = flip(LaHs1RHO,3);   
  vvel_tavg = VVEL;
  temp_tavg = THETA;
  salt_tavg = SALT;  
  pressure_tavg =PHIHYD; 

  temp_tavg(hFacC==0) = NaN;
  salt_tavg(hFacC==0) = NaN;
  %%% Calculate the potential density pt
  refdepth = -zz(layers_krho(1)); %%% equals: sum(DRF(1:layers_krho-1))+DRF(layers_krho)/2;
  rhoConst = 1027.0 %??? is this correct?
  g=9.8;
  refpress = rhoConst*(g*refdepth + pressure_tavg(:,:,layers_krho(1)))/1e4; %%% unit: dbar
%   refpress = refdepth;
  
  for kk = 1:Nr
      pt_tavg(:,:,kk) = densmdjwf(salt_tavg(:,:,kk),temp_tavg(:,:,kk),refpress)-1000;
  end
  
  pt_tavg(hFacC==0) = NaN;

%   vflux_tavg = zeros(Nx,Ny,Npt);
%   h_pt_tavg = zeros(Nx,Ny,Npt);
%   pt_tavg = zeros(Nx,Ny,Nr);
%   vvel_tavg = zeros(Nx,Ny,Nr);
%   navg = 0;
%   for n=1:length(dumpIters)
% 
%     tyears = dumpIters(n)*deltaT/86400/365;
% 
%     if ((tyears >= tmin) && (tyears <= tmax))    
% 
%       [tyears dumpIters(n)]
%       vflux = rdmdsWrapper(fullfile(exppath,'results/LaVH1TH'),dumpIters(n));      
%       h_pt = rdmdsWrapper(fullfile(exppath,'results/LaHs1TH'),dumpIters(n));      
%       pt  = rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters(n));                      
%       vvel  = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),dumpIters(n));              
% 
%       if (isempty(vflux) || isempty(pt) ...
%           || isempty(h_pt) || isempty(vvel))
%         ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tyears),' days.']
%         break;
%       else
%         vflux_tavg = vflux_tavg + vflux;
%         h_pt_tavg = h_pt_tavg + h_pt;
%         pt_tavg = pt_tavg + squeeze(pt(:,:,:,1));      
%         vvel_tavg = vvel_tavg + squeeze(vvel(:,:,:,1));      
%         navg = navg + 1;
%       end
%     end
% 
%   end
% 
%   %%% Calculate the time average
%   if (navg == 0)
%     error('No data files found');
%   end
%   vflux_tavg = vflux_tavg/navg;
%   h_pt_tavg = h_pt_tavg/navg;
%   pt_tavg = pt_tavg/navg;
%   vvel_tavg = vvel_tavg/navg;
%   pt_tavg(hFacC==0) = NaN;
 
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
  flip_ptlevs = flip(ptlevs);
  vflux_m(:,:,Npt) = vflux_m(:,:,Npt) + sum(vdz.*(pt_f>flip_ptlevs(Npt)),3);
  vflux_m(:,:,1) = vflux_m(:,:,1) + sum(vdz.*(pt_f<=flip_ptlevs(2)),3);
  for m=2:Npt-1
    vflux_m(:,:,m) = vflux_m(:,:,m) + sum(vdz.*((pt_f>flip_ptlevs(m)) & (pt_f<=flip_ptlevs(m+1))),3);
  end   
 
  %%% Zonally integrate meridional fluxes
  vflux_xint = zeros(Ny,Npt);
  vflux_m_xint = zeros(Ny,Npt);
  for i=1:Nx
    vflux_xint = vflux_xint + delX(i)*squeeze(vflux_tavg(i,:,:));
    vflux_m_xint = vflux_m_xint + delX(i)*squeeze(vflux_m(i,:,:));
  end
 
  vflux_m_xint = flip(vflux_m_xint,2);

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
    z_pt(:,m) = - sum(h_pt_xtavg(:,1:m-1),2);
  end
 
  %%% Calculate zonal-mean potential temperature
  pt_xtavg = squeeze(nanmean(pt_tavg(:,:,:)));
  pt_f_xtavg = squeeze(nanmean(pt_f(:,:,:)));
 
  %%% Convert to z-coordinates by mapping the streamfunction at each temp 
  %%% level to the mean height of that density surface
  psi_z = NaN*ones(Ny,Nrf);
  psim_z = NaN*ones(Ny,Nrf);
  psie_z = NaN*ones(Ny,Nrf);
  for j=1:Ny  
 
    for k=1:Nrf
 
      %%% Density lies in the lowest bin
      if (pt_f_xtavg(j,k) > ptlevs(2))
        psi_z(j,k) = psi_pt(j,1);      
        psim_z(j,k) = psim_pt(j,1);    
        psie_z(j,k) = psie_pt(j,1);    
        continue;
      end
 
      %%% Density lies in the highest bin
      if (pt_f_xtavg(j,k) < ptlevs(Npt))
        psi_z(j,k) = psi_pt(j,Npt);      
        psim_z(j,k) = psim_pt(j,Npt);      
        psie_z(j,k) = psie_pt(j,Npt);      
        continue;
      end    
 
      %%% Density lies in an intermediate bin, so find the bin and assign
      %%% the overturning streamfunction via linear interpolation
      for m=2:Npt-1
        if (pt_f_xtavg(j,k) > ptlevs(m+1))
          pt_n = ptlevs(m+1);
          pt_p = ptlevs(m);
          wp = (pt_n-pt_f_xtavg(j,k))/(pt_n-pt_p);
          wn = 1 - wp;
          psi_z(j,k) = wp*psi_pt(j,m) + wn*psi_pt(j,m+1);
          psim_z(j,k) = wp*psim_pt(j,m) + wn*psim_pt(j,m+1);
          psie_z(j,k) = wp*psie_pt(j,m) + wn*psie_pt(j,m+1);
          break;
        end
      end
 
    end
 
  end

  

 
  figure(10)
    PSIlim=[-1.5 1.5];    YLIM = [36.4 37.5];
    subplot(1,3,1)
    pcolor(yy/1000,ptlevs,psi_pt');
    shading interp;colormap('redblue');colorbar;caxis(PSIlim);axis ij;
    ylim(YLIM);
    xlabel('y (km)');ylabel('Potential density (kg/m^3)');
    title('\psi (Sv) in PT space')
    subplot(1,3,2)
    pcolor(yy/1000,ptlevs,psim_pt');
    shading interp;colormap('redblue');colorbar;caxis(PSIlim);axis ij;
    xlabel('y (km)');ylabel('z (m)');    ylim(YLIM);
    title('\psi_{mean} (Sv) in PT space')
    subplot(1,3,3)
    pcolor(yy/1000,ptlevs,psie_pt');
    shading interp;colormap('redblue');colorbar;caxis(PSIlim);axis ij;
    xlabel('y (km)');ylabel('z (m)');    ylim(YLIM);
    title('\psi_{eddy} (Sv) in PT space')
    
  figure(11);
    PSIlim=[-1.5 1.5];    YLIM = [36.4 37.5];
    subplot(1,3,1)
    pcolor(yy/1000,zz_f,psi_z');
    shading interp;colormap('redblue');colorbar;caxis(PSIlim);
  
    xlabel('y (km)');ylabel('Potential density (kg/m^3)');
    title('\psi (Sv) in z space')
    subplot(1,3,2)
    pcolor(yy/1000,zz_f,psim_z');
    shading interp;colormap('redblue');colorbar;caxis(PSIlim);
    xlabel('y (km)');ylabel('z (m)');   
    title('\psi_{mean} (Sv) in z space')
    subplot(1,3,3)
    pcolor(yy/1000,zz_f,psie_z');
    shading interp;colormap('redblue');colorbar;caxis(PSIlim);
    xlabel('y (km)');ylabel('z (m)');    
    title('\psi_{eddy} (Sv) in z space')
  %%

  %%% Store computed data for later
  save(fullfile(prodir,[expname,'_MOC_rho.mat']), ...
    'xx','yy','zz','zz_f','ptlevs', ... 
    'pt_tavg','pt_xtavg','pt_f_xtavg', ...  
    'vflux_m','vflux_xint','vflux_m_xint', ...
    'psi_pt','psim_pt','psie_pt', ...
    'psi_z','psim_z','psie_z');
  
  
end
