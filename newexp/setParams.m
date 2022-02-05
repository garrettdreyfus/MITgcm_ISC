%%%
%%% setParams.m
%%%
%%% Sets basic MITgcm parameters plus parameters for included packages, and
%%% writes out the appropriate input files.,
%%%
function nTimeSteps = setParams (exp_name,inputpath,codepath,listterm,Nx,Ny,Nr)

  %%% Load EOS utilities
  addpath ~/UCLA/Utilities/GSW/
  addpath ~/UCLA/Utilities/GSW/library/
  
  %%% Other required tools
  addpath ../utils/matlab;
  addpath ../newexp_utils/;
  
  
  
  
  %%%%%%%%%%%%%%%%%%
  %%%%% SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%      
  
  %%% If set true, plots of prescribed variables will be shown
  showplots = true;      
  fignum = 1;
  fontsize = 12;
  
  %%% Data format parameters
  ieee='b';
  prec='real*8';
  realdigits = 8;
  realfmt=['%.',num2str(realdigits),'e'];
  
  %%% Get parameter type definitions
  paramTypes;
      

  %%% To store parameter names and values
  parm01 = parmlist;
  parm02 = parmlist;
  parm03 = parmlist;
  parm04 = parmlist;
  parm05 = parmlist;
  PARM={parm01,parm02,parm03,parm04,parm05}; 
  
  %%% Seconds in one hour
  t1min = 60;
  %%% Seconds in one hour
  t1hour = 60*t1min;
  %%% hours in one day
  t1day = 24*t1hour;
  %%% Seconds in 1 year
  t1year = 365*t1day;  
  %%% Metres in one kilometre
  m1km = 1000; 
      
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% FIXED PARAMETER VALUES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  simTime = 20*t1year; %%% Simulation time   
%   simTime = 60*t1day;
  nIter0 = 0; %%% Initial iteration 
  Lx = 400*m1km; %%% Domain size in x 
  Ly = 400*m1km; %%% Domain size in y   
%   Ls = 50*m1km; %%% Width of southern boundary region
  Ln = 20*m1km; %%% Width of northern boundary region
  H = 4000; %%% Domain size in z 
  g = 9.81; %%% Gravity  
  f0 = -1.3e-4; %%% Coriolis parameter
  beta = 1e-11; %%% Beta parameter      
  rhoConst = 1027; %%% Reference density
  
  viscAh = 0; %%% Horizontal viscosity    
  viscA4 = 0; %%% Biharmonic viscosity
  viscAhGrid = 0; %%% Grid-dependent viscosity
  viscAr = 1e-4; %%% Vertical viscosity    
  diffKhT = 0; %%% Horizontal temp diffusion
  diffKhS = 0; %%% Horizontal salt diffusion
%   viscA4Grid = 0.1; %%% Grid-dependent biharmonic viscosity
%   viscC4smag = 0; %%% Smagorinsky biharmonic viscosity  
%   diffK4Tgrid = 0.1; %%% Grid-dependent biharmonic temp diffusivity
  diffKrT = 1e-5; %%% Vertical temp diffusion     
%   diffK4Sgrid = 0.1; %%% Grid-dependent biharmonic salt diffusivity
  diffKrS = 1e-5; %%% Vertical salt diffusion     
  viscA4Grid = 0;    %%%%% Update: 20210627
  viscC4smag = 4;    %%%%% Update: 20210627
  diffK4Tgrid = 0;   %%%%% Update: 20210627
%   diffKrT = 0;       %%%%% Update: 20210702
  diffK4Sgrid = 0;   %%%%% Update: 20210627
%   diffKrS = 0;       %%%%% Update: 20210702
%   viscAh = 12;       %%%%% Update: 20210630
%   viscA4Grid = 0.1;  %%%%% Update: 20210630

  
  
  %%% Topographic parameters 
  Wslope = 30*m1km; %%% Continental slope half-width
  Hshelf = 500; %%% Continental shelf depth
  Wshelf = 50*m1km; %%% Width of continental shelf
  Ycoast = 200*m1km; %%% Latitude of coastline
  Wcoast = 20*m1km; %%% Width of coastal wall slope
  Yshelfbreak = Ycoast+Wshelf; %%% Latitude of shelf break
  Yslope = Ycoast+Wshelf+Wslope; %%% Latitude of mid-continental slope
  Ydeep = 350*m1km; %%% Latitude of deep ocean
  Xeast = 275*m1km; %%% Longitude of eastern trough wall
  Xwest = 125*m1km; %%% Longitude of western trough wall
  Yicefront = 150*m1km; %%% Latitude of ice shelf face
  Hicefront = 200; %%% Depth of ice shelf frace
  Hbed = -500; %%% Change in bed elevation from shelf break to southern domain edge
  Hice = Hicefront-(Hshelf-Hbed); %%% Change in ice thickness from ice fromt to southern domain edge
  Htrough = 300; %%% Trough depth
  Wtrough = 15*m1km; %%% Trough width
  Xtrough = Lx/2; %%% Longitude of trough

  %%% Package options
  useSEAICE = false;
  useSHELFICE = true;     
  useLAYERS = false;  
  
  %%% TODO  
  useEXF = useSEAICE;
  EXFoption = 3; %%% Read-in atemp, aqh, swdown, lwdown,precip, and runoff. Compute hflux, swflux and sflux. 
  
  %%% OBCS package options
  useOBCS = true;    
  useOBCSbalance = true;  
  useOrlanskiNorth = false;
  
  
  %%% PARM01
  %%% momentum scheme
  parm01.addParm('vectorInvariantMomentum',true,PARM_BOOL);
  parm01.addParm('implicSurfPress',0.6,PARM_REAL);
  parm01.addParm('implicDiv2DFlow',0.6,PARM_REAL); 
  %%% viscosity  
  parm01.addParm('viscAr',viscAr,PARM_REAL);
  parm01.addParm('viscA4',viscA4,PARM_REAL);
  parm01.addParm('viscAh',viscAh,PARM_REAL);
  parm01.addParm('viscA4Grid',viscA4Grid,PARM_REAL);
  parm01.addParm('viscAhGrid',viscAhGrid,PARM_REAL);
  parm01.addParm('viscA4GridMax',0.5,PARM_REAL);
  parm01.addParm('viscAhGridMax',1,PARM_REAL);
  parm01.addParm('useAreaViscLength',false,PARM_BOOL);
  parm01.addParm('useFullLeith',true,PARM_BOOL);
  parm01.addParm('viscC4smag',viscC4smag,PARM_REAL);    
  parm01.addParm('viscC4leith',0,PARM_REAL);
  parm01.addParm('viscC4leithD',0,PARM_REAL);  
  parm01.addParm('viscC2leith',0,PARM_REAL);
  parm01.addParm('viscC2leithD',0,PARM_REAL);  
  %%% diffusivity
  parm01.addParm('tempAdvScheme',80,PARM_INT);
  parm01.addParm('saltAdvScheme',80,PARM_INT);
  parm01.addParm('diffKrT',diffKrT,PARM_REAL);
  parm01.addParm('diffKhT',diffKhT,PARM_REAL);
  parm01.addParm('diffKrS',diffKrS,PARM_REAL);
  parm01.addParm('diffKhS',diffKhS,PARM_REAL);
  parm01.addParm('tempStepping',true,PARM_BOOL);
  parm01.addParm('saltStepping',true,PARM_BOOL);
  parm01.addParm('staggerTimeStep',true,PARM_BOOL);
  %%% equation of state
  parm01.addParm('eosType','MDJWF',PARM_STR); 
  %%% boundary conditions
  parm01.addParm('no_slip_sides',false,PARM_BOOL);
  parm01.addParm('no_slip_bottom',false,PARM_BOOL);
  parm01.addParm('bottomDragLinear',0,PARM_REAL);
  parm01.addParm('bottomDragQuadratic',2e-3,PARM_REAL);
  %%% physical parameters
  parm01.addParm('f0',f0,PARM_REAL);
  parm01.addParm('beta',beta,PARM_REAL);
  parm01.addParm('gravity',g,PARM_REAL);
  %%% full Coriolis force parameters
  parm01.addParm('quasiHydrostatic',false,PARM_BOOL);
  parm01.addParm('fPrime',0,PARM_REAL);
  parm01.addParm('rhoConst',rhoConst,PARM_REAL);
  %%% implicit diffusion and convective adjustment  
  parm01.addParm('ivdc_kappa',0,PARM_REAL);
  parm01.addParm('implicitDiffusion',true,PARM_BOOL);
  parm01.addParm('implicitViscosity',true,PARM_BOOL);
  %%% exact volume conservation
  parm01.addParm('exactConserv',true,PARM_BOOL);
  %%% C-V scheme for Coriolis term
  parm01.addParm('useCDscheme',false,PARM_BOOL);
  %%% partial cells for smooth topography  
  parm01.addParm('hFacMin',0.1,PARM_REAL);    
  %%% file IO stuff
  parm01.addParm('readBinaryPrec',64,PARM_INT);
  parm01.addParm('useSingleCpuIO',true,PARM_BOOL);
  parm01.addParm('debugLevel',2,PARM_INT);
  %%% Wet-point method at boundaries - may improve boundary stability
  parm01.addParm('useJamartWetPoints',true,PARM_BOOL);
  parm01.addParm('useJamartMomAdv',true,PARM_BOOL);
%   parm01.addParm('rhoConst',1000,PARM_REAL);
  parm01.addParm('useRealFreshWaterFlux',false,PARM_BOOL);
  %%% PARM02
  parm02.addParm('useSRCGSolver',true,PARM_BOOL);  
  parm02.addParm('cg2dMaxIters',1000,PARM_INT);  
  parm02.addParm('cg2dTargetResidual',1e-12,PARM_REAL);
 
  %%% PARM03
  parm03.addParm('alph_AB',1/2,PARM_REAL);
  parm03.addParm('beta_AB',5/12,PARM_REAL);
  if (useSEAICE)
      parm03.addParm('forcing_In_AB',false,PARM_BOOL); 
      % This flag makes to model do a  separate (Eulerian?) time step 
      % for the tendencies due to surface forcing. This is sometime 
      % favorable for stability reasons (and some package such as 
      % seaice work only with this).
  end
  parm03.addParm('nIter0',nIter0,PARM_INT);
  parm03.addParm('abEps',0.1,PARM_REAL);
  parm03.addParm('chkptFreq',t1year/12,PARM_REAL); % rolling 
  parm03.addParm('pChkptFreq',t1year,PARM_REAL); % permanent
  parm03.addParm('taveFreq',0,PARM_REAL); % it only works properly, if taveFreq is a multiple of the time step deltaT (or deltaTclock).
  parm03.addParm('dumpFreq',0,PARM_REAL); % interval to write model state/snapshot data (s)
  parm03.addParm('monitorFreq',t1year/12,PARM_REAL); % interval to write monitor output (s)
  parm03.addParm('dumpInitAndLast',true,PARM_BOOL);
  parm03.addParm('pickupStrictlyMatch',false,PARM_BOOL); 
  
  %%% PARM04
  parm04.addParm('usingCartesianGrid',true,PARM_BOOL);
%   parm04.addParm('usingCurvilinearGrid',true,PARM_BOOL);
  parm04.addParm('usingSphericalPolarGrid',false,PARM_BOOL);    
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% GRID SPACING %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%    

  %%% Zonal grid
  dx = Lx/Nx;  
  xx = (1:Nx)*dx;
  
  %%% Uniform meridional grid   
  dy = (Ly/Ny)*ones(1,Ny);  
  yy = cumsum((dy + [0 dy(1:end-1)])/2);
 
  %%% Plotting mesh
  [Y,X] = meshgrid(yy,xx);
   
%   %%% Grid spacing increases with depth, but spacings exactly sum to H
%   zidx = 1:Nr;
%   gamma = 10;  
%   alpha = 100;
%   dz1 = 2*H/Nr/(alpha+1);
%   dz2 = alpha*dz1;
%   dz = dz1 + ((dz2-dz1)/2)*(1+tanh((zidx-((Nr+1)/2))/gamma));
%   zz = -cumsum((dz+[0 dz(1:end-1)])/2);

  %%% Variable grid with high resolution at ice shelf cavity depths, very high in surface mixed layer    
  dz0 = 2;
  dz1 = 15; 
  dz2 = 20;
  dz3 = 100;
  dz4 = 200;  
  N0 = 1;
  N1 = 20; 
  N2 = 50;
  N3 = 15;
  N4 = 14;  
  nn_c = cumsum([N0 N1 N2 N3 N4]);
  dz_c = [dz0 dz1 dz2 dz3 dz4];
  nn = 1:(N1+N2+N3+N4+1);
  dz = interp1(nn_c,dz_c,nn,'pchip');

  zz = -cumsum((dz+[0 dz(1:end-1)])/2);
  if (length(zz) ~= Nr)
    error('Vertical grid size does not match vertical array dimension!');
  end
  Nr = length(zz);


  %%% Thickness of sponge layers in gridpoints  
  spongeThicknessDim = Ln;
  spongeThickness = round(spongeThicknessDim/dy(end));
  seaicespongeThicknessDim = Ln; %%% Restore sea ice thickness and concentration for the whole domain
  seaiceSpongeThickness = round(seaicespongeThicknessDim/dy(end));
  

  %%% Store grid spacings
  parm04.addParm('delX',dx*ones(1,Nx),PARM_REALS);
  parm04.addParm('delY',dy,PARM_REALS);
  parm04.addParm('delR',dz,PARM_REALS);      
  
  %%% Don't allow partial cell height to fall below min grid spacing  
  parm01.addParm('hFacMinDr',min(dz),PARM_REAL);  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% BATHYMETRY AND ICE SHELF DRAFT %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  %%% Construct shelf/slope/deep ocean bathymetry profile via cubic
  %%% interpolation
  y_interp = [0 (Yslope-Wslope)/2 Yslope-Wslope Yslope Ydeep Ly];
  h_interp = [-Hshelf+Hbed -Hshelf+Hbed/2 -Hshelf -(Hshelf+H)/2 -H -H];
  h_profile = interp1(y_interp,h_interp,yy,'pchip');
  h = repmat(h_profile,[Nx 1]);
  
  
  %%% Add trough
  y_interp = [0 Yshelfbreak Yslope Ly];
  h_interp = [0 -Htrough 0 0];
  h_trough_profile = interp1(y_interp,h_interp,yy,'pchip');
  h_trough = repmat(h_trough_profile,[Nx 1]);
  h_trough = h_trough .* 1./(cosh((X-Xtrough)/Wtrough)).^2;
  h = h + h_trough;
  
    
  %%% Add coastal wall %%%
  h_coast = zeros(Nx,Ny);
  
  %%% Western coastline
  coastidx = (Y<Ycoast+Wcoast/2) & (Y>Ycoast-Wcoast/2) & (X<=Xwest-Wcoast/2);
  h_coast(coastidx) = h_coast(coastidx) - h(coastidx).*0.5.*(1-sin(pi*(Y(coastidx)-Ycoast)/Wcoast));
  landidx = find((Y<=Ycoast-Wcoast/2) & (X<=Xwest-Wcoast/2));
  h_coast(landidx) = -h(landidx);
  
  %%% Western corner
  R = sqrt((X-(Xwest-Wcoast/2)).^2+(Y-(Ycoast-Wcoast/2)).^2);
  coastidx = (Y>Ycoast-Wcoast/2) & (X>Xwest-Wcoast/2) & (R <= Wcoast);  
  h_coast(coastidx) = h_coast(coastidx) - h(coastidx).*0.5.*(1+cos(pi*R(coastidx)/Wcoast));
  
  %%% Western trough wall
  coastidx = (Y<Ycoast-Wcoast/2) & (X<=Xwest+Wcoast/2) & (X>Xwest-Wcoast/2);
  h_coast(coastidx) = h_coast(coastidx) - h(coastidx).*0.5.*(1-sin(pi*(X(coastidx)-Xwest)/Wcoast));   
  
  %%% Eastern coastline
  coastidx = (Y<Ycoast+Wcoast/2) & (Y>Ycoast-Wcoast/2) & (X>=Xeast+Wcoast/2);
  h_coast(coastidx) = h_coast(coastidx) - h(coastidx).*0.5.*(1-sin(pi*(Y(coastidx)-Ycoast)/Wcoast));
  landidx = find((Y<=Ycoast-Wcoast/2) & (X>=Xeast+Wcoast/2));
  h_coast(landidx) = -h(landidx);
  
  %%% Eastern corner
  R = sqrt((X-(Xeast+Wcoast/2)).^2+(Y-(Ycoast-Wcoast/2)).^2);
  coastidx = (Y>Ycoast-Wcoast/2) & (X<Xeast+Wcoast/2) & (R <= Wcoast);  
  h_coast(coastidx) = h_coast(coastidx) - h(coastidx).*0.5.*(1+cos(pi*R(coastidx)/Wcoast));
  
  %%% Eastern trough wall
  coastidx = (Y<Ycoast-Wcoast/2) & (X<Xeast+Wcoast/2) & (X>=Xeast-Wcoast/2);
  h_coast(coastidx) = h_coast(coastidx) - h(coastidx).*0.5.*(1+sin(pi*(X(coastidx)-Xeast)/Wcoast));   
  
  h = h + h_coast;
 
  
  
 
  %%% Construct ice shelf
  icedraft = zeros(Nx,Ny);
  iceidx = find(Y<=Yicefront);  
  icedraft(iceidx) = -Hicefront - (Y(iceidx)-Yicefront)/Yicefront * Hice;
  icedraft(icedraft<h) = h(icedraft<h);
  
  
  %%% Make sure there are no "holes" along the southern boundary, or MITgcm
  %%% will think it's supposed to be north/south periodic
  wallidx = find(icedraft(:,1)>h(:,1));
  h(wallidx,1) = icedraft(wallidx,1);
  
  
  %%% Plot bathymetry and ice draft
  figure(fignum);
  fignum = fignum + 1;
  clf;    

  %%% Bathymetry  
  p = surface(X(:,2:end-1)/1000,Y(:,2:end-1)/1000,h(:,2:end-1));
  p.FaceColor = [11*16+9 9*16+12 6*16+11]/255;
  p.EdgeColor = 'none';

  %%% Modified ice draft to look good in the plot
  icedraft_plot = icedraft;
  icedraft_plot(icedraft==0) = NaN;
  icetop_plot = 0*icedraft_plot;
  for i=1:Nx
    j = find(~isnan(icetop_plot(i,:)),1,'last');
    if (isempty(j))
      continue;
    else
      icetop_plot(i,j+1) = max(-Hicefront,h(i,j+1));
    end
  end
 
  %%% Plot ice
  hold on;
  p = surface(X(:,2:end-1)/1000,Y(:,2:end-1)/1000,icedraft_plot(:,2:end-1));
  p.FaceColor = [153, 255, 255]/255;
  p.EdgeColor = 'none';
  p = surface(X(:,2:end-1)/1000,Y(:,2:end-1)/1000,icetop_plot(:,2:end-1));
  p.FaceColor = [153, 255, 255]/255;
  p.EdgeColor = 'none';
  hold off;

  %%% Decorations
  view(-206,14);
  axis tight;
  xlabel('x (km)','interpreter','latex');
  ylabel('y (km)','interpreter','latex');
  zlabel('z (m)','interpreter','latex');
  pbaspect([Lx/Ly 1 1]);
  camlight('headlight');
  lightangle(-206,34);
  lighting gouraud;
  
  %%% Save as a parameter
  writeDataset(h,fullfile(inputpath,'bathyFile.bin'),ieee,prec);
  parm05.addParm('bathyFile','bathyFile.bin',PARM_STR); 
  
  
  
 
  
  
  %%% TODO buoyancy frequency looks very low close to the surface
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% NORTHERN TEMPERATURE/SALINITY PROFILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Bottom properties offshore, taken from Meijers et al. (2010)
  %%% measurements. We need these because the KN climatology only goes down
  %%% to 2000m
  s_bot = 34.66;
  pt_bot = -0.5;
  
%   %%% Load sections from Kapp Norvegia climatology (Hattermann 2018) 
  ptemp_KN = ncread('KappNorvegiaCLM.nc','ptemp'); %%% at sea level pressure
  salt_KN = ncread('KappNorvegiaCLM.nc','salt');
  pres_KN = ncread('KappNorvegiaCLM.nc','pressure');
    
  
  %%% Winter means at offshore boundary
  ptemp_North = [squeeze(mean(ptemp_KN(:,end,6:8),3))' pt_bot];
  salt_North = [squeeze(mean(salt_KN(:,end,6:8),3))' s_bot];
%   N_offshore = 40:45;
%   ptemp_North = [squeeze(mean(mean(ptemp_KN(:,N_offshore,6:8),3),2))' pt_bot];
%   salt_North = [squeeze(mean(mean(salt_KN(:,N_offshore,6:8),3),2))' s_bot];
  depth_North = [-pres_KN' -H];
  
  %%% Interpolate to model grid
  tNorth = interp1(depth_North,ptemp_North,zz,'PCHIP'); %%% reference pressure level: sea surface
  sNorth = interp1(depth_North,salt_North,zz,'PCHIP');  %%% reference pressure level: sea surface
  
  
%   zidx = find(sNorth==max(sNorth));
%   %   sNorth(zidx+1:end) = 0.5*(sNorth(zidx+1:end)+max(sNorth));
%   sNorth(zidx+1:end) = max(sNorth);
%   tNorth(1) = -1.87;
  
%   load('NorthernBdry')
%   sNorth = Snorth_exp_dense';
%   tNorth = Tnorth_exp_dense';
%   tNorth(1) = -1.87;

%%% Calculate the freezing temperature at ocean surface using GSW toolbox
%   Ai0 = 1; % Initial fractional sea ice cover,
%   Hi0 = 1; % Initial sea ice thickness = 1m
  Hs0 = 0; % Initial snow thickness = 0.1 m
  Si0 = 6; % The salinity for 1m sea ice is about 6 g/kg. Cox et al., (1974). Salinity variations in sea ice. 
  rho_i = 920; % Density of sea ice
  pp0 = - zz; % This pressure is approximate, using a constant density

  pp = pp0;
  % %%% Calculate the true hydrostatic pressure, need to check
  %     pp = zeros(size(zz));
  %     DZ = zeros(size(zz));
  %     DZ(1) = zz(1);
  %     DZ(2:end) = zz(2:end)-zz(1:end-1);
  %     insitu_rho_north_densmdjwf = densmdjwf(sNorth,tNorth,pp0)';         
  %     P0 = 0; % surface atmospheric pressure
  %     pp(1) = P0 - insitu_rho_north_densmdjwf(1)*g*DZ(1)/10^4;
  %     for nz = 2:size(zz,2)
  %         insitu_rho_north_densmdjwf(nz) = densmdjwf(sNorth(nz),tNorth(nz),pp(nz-1));
  %         pp(nz) = pp(nz-1) - insitu_rho_north_densmdjwf(nz)*g*DZ(nz)/10^4;
  %     end
    
  %%% Calculate the freezing temperature using GSW toolbox
  %   pp = pp0;
  %   pfreezing = rho_i*g*Hi0*Ai0/10^4+pp(1);
  %   SA0 = gsw_SA_from_SP(sNorth(1),pfreezing,-12,-64);
  %   saturation_fraction = 1 - Ai0; 
  %   tNorth(1) = gsw_t_freezing(SA0,pfreezing,saturation_fraction);  %%% -1.8862
  %   end 


    
  
  %%% Plot the relaxation temperature
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(tNorth,-zz,'LineWidth',1.5); axis ij;    
    xlabel('\theta_r_e_f (\circC)');
    ylabel('Depth (m)');
    title('Relaxation temperature');
    legend('Northern T','Southern T','Position',[0.3200 0.6468 0.3066 0.0738]);
    set(gca,'fontsize',fontsize);
    PLOT = gcf;
    PLOT.Position = [644 148 380 562];  
  end
    
  %%% Plot the relaxation salinity
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(sNorth,-zz,'LineWidth',1.5);axis ij;
    xlabel('S_r_e_f (psu)');
    ylabel('Depth (m)');
%     ylabel('z','Rotation',0);
    title('Relaxation salinity');
    legend('Northern S','Southern S','Position',[0.3200 0.6468 0.3066 0.0738]);
    set(gca,'fontsize',fontsize);
    PLOT = gcf;
    PLOT.Position = [644 148 380 562];  
  end

  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DEFORMATION RADIUS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
  %%% Check Brunt-Vaisala frequency using full EOS
  %%% TODO need to calculate SA and CT
  SA_north = sNorth;
  CT_north = tNorth;
  [N2_north, pp_mid_north] = gsw_Nsquared(SA_north,CT_north,pp,-64);  
  dzData = zz(1:end-1)-zz(2:end);

  %%% Calculate internal wave speed and first Rossby radius of deformation
  N = sqrt(N2_north);
  Cig = zeros(size(yy));
  for j=1:Ny    
    for k=1:length(dzData)
      if (zz(k) > h(1,j))        
        Cig(j) = Cig(j) + N(k)*min(dzData(k),zz(k)-h(1,j));
      end
    end
  end
  Rd = Cig./(pi*abs(f0+beta*Y(1,:)));

  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    semilogx(N2_north,pp_mid_north,'LineWidth',1.5);axis ij;
    legend('Northern N^2','Southern N^2','Position',[0.5181 0.6192 0.3313 0.0899]);
    xlabel('N^2 (s^-^2)');
    ylabel('Depth (m)');
%       ylabel('z (km)','Rotation',0);
    title('Buoyancy frequency');
    set(gca,'fontsize',fontsize);
    PLOT = gcf;
    PLOT.Position = [644 148 380 562];  
  end

  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(yy/1000,Rd/1000,'LineWidth',1.5);
    xlabel('Offshore distance (km)');
    ylabel('R_d (km)');
    title('First baroclinic Rossby deformation radius');
    set(gca,'fontsize',fontsize-1);
    PLOT = gcf;
    PLOT.Position = [263 149 567 336];  
  end  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CALCULATE TIME STEP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
  
  %%% These estimates are in no way complete, but they give at least some
  %%% idea of the time step needed to keep things stable. In complicated 
  %%% simulations, preliminary tests may be required to estimate the
  %%% parameters used to calculate these time steps.        
  
  %%% Gravity wave CFL

  %%% Upper bound for absolute horizontal fluid velocity (m/s)
  %%% At the moment this is just an estimate
%   Umax = 1
  Umax = 2;
  %%% Max gravity wave speed
  cmax = max(Cig);
  %%% Max gravity wave speed using total ocean depth
  cgmax = Umax + cmax;
  %%% Advective CFL
  deltaT_adv = min([0.5*dx/cmax,0.5*dy/cmax]);
  %%% Gravity wave CFL
  deltaT_gw = min([0.5*dx/Umax,0.5*dy/Umax]);
  %%% CFL time step based on full gravity wave speed
  deltaT_fgw = min([0.5*dx/cgmax,0.5*dy/cgmax]);
    
  %%% Other stability conditions
  
  %%% Inertial CFL time step (Sf0<=0.5)
  deltaT_itl = 0.5/abs(f0);
  %%% Time step constraint based on horizontal diffusion 
  deltaT_Ah = 0.5*min([dx dy])^2/(4*viscAh);    
  %%% Time step constraint based on vertical diffusion
  deltaT_Ar = 0.5*min(dz)^2 / (4*viscAr);  
  %%% Time step constraint based on biharmonic viscosity 
  deltaT_A4 = 0.5*min([dx dy])^4/(32*viscA4);
  %%% Time step constraint based on horizontal diffusion of temp 
  deltaT_KhT = 0.4*min([dx dy])^2/(4*diffKhT);    
  %%% Time step constraint based on vertical diffusion of temp 
  deltaT_KrT = 0.4*min(dz)^2 / (4*diffKrT);
  
  %%% Time step size  
  deltaT = min([deltaT_fgw deltaT_gw deltaT_adv deltaT_itl deltaT_Ah deltaT_Ar deltaT_KhT deltaT_KrT deltaT_A4]);
  deltaT = round(deltaT) 
%   deltaT = 100
  nTimeSteps = ceil(simTime/deltaT);
  simTimeAct = nTimeSteps*deltaT
  
  %%% Write end time time step size  
  parm03.addParm('endTime',nIter0*deltaT+simTimeAct,PARM_INT);
  parm03.addParm('deltaT',deltaT,PARM_REAL); 

  
 
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL DATA %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%
    
  %%% Random noise amplitude
  tNoise = 0.001;  
  sNoise = 0.0001;
      
  %%% Align initial temp with background
  hydroTh = ones(Nx,Ny,Nr);
  hydroSa = ones(Nx,Ny,Nr);
  for k=1:1:Nr
    hydroTh(:,:,k) = squeeze(hydroTh(:,:,k))*tNorth(k);
    hydroSa(:,:,k) = squeeze(hydroSa(:,:,k))*sNorth(k);
  end
  
  %%% Add some random noise  
  hydroTh = hydroTh + tNoise*(2*rand(Nx,Ny,Nr)-1);
  hydroSa = hydroSa + sNoise*(2*rand(Nx,Ny,Nr)-1);  
  
  %%% Write to data files
  writeDataset(hydroTh,fullfile(inputpath,'hydrogThetaFile.bin'),ieee,prec); 
  parm05.addParm('hydrogThetaFile','hydrogThetaFile.bin',PARM_STR);
  writeDataset(hydroSa,fullfile(inputpath,'hydrogSaltFile.bin'),ieee,prec); 
  parm05.addParm('hydrogSaltFile','hydrogSaltFile.bin',PARM_STR); 
  
  
    
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% TRACER DIFFUSION %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Set biharmonic diffusivities as a fraction of the maximum stable
  %%% grid-scale hyperdiffusion
  diffK4T = diffK4Tgrid * max([dx dy])^4 / (32*deltaT);
  diffK4S = diffK4Sgrid * max([dx dy])^4 / (32*deltaT);
  parm01.addParm('diffK4T',diffK4T,PARM_REAL); 
  parm01.addParm('diffK4S',diffK4S,PARM_REAL); 
  
  
  
  
  
 
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates the 'data' file
  write_data(inputpath,PARM,listterm,realfmt);
  
 

  
  
  
  
  
  
  

  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%% SHELF ICE   %%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  if (useSHELFICE)
 
    %%% to store parameter names and values
    shelfice_parm01 = parmlist;
    SHELFICE_PARM = {shelfice_parm01}; 

    SHELFICEHeatCapacity_Cp = 2000; %%% Default value
    rhoShelfIce = 917;              %%% Default value
    SHELFICEheatTransCoeff = 0;     %%% Turn off linear heat transfer
    SHELFICEthetaSurface = -20;     %%% Default value %%% TODO should match atmospheric temp over the sea ice?
    SHELFICEuseGammaFrict = true;   %%% Turn on friction-dependent heat transfer
    SHELFICEboundaryLayer = false;  %%% Turn on to average velocities over top dz of water column when computing friction velocity
    SHELFICEconserve = false;       %%% Turns on conservative form of 3-equation IOBL parameterization
    
    %%% Save as a parameter   
    shelfice_parm01.addParm('SHELFICEHeatCapacity_Cp',SHELFICEHeatCapacity_Cp,PARM_REAL);
    shelfice_parm01.addParm('rhoShelfIce',rhoShelfIce,PARM_REAL);
    shelfice_parm01.addParm('SHELFICEheatTransCoeff',SHELFICEheatTransCoeff,PARM_REAL);
    shelfice_parm01.addParm('SHELFICEthetaSurface',SHELFICEthetaSurface,PARM_REAL);
    shelfice_parm01.addParm('SHELFICEuseGammaFrict',SHELFICEuseGammaFrict,PARM_BOOL);
    shelfice_parm01.addParm('SHELFICEboundaryLayer',SHELFICEboundaryLayer,PARM_BOOL);
    shelfice_parm01.addParm('SHELFICEconserve',SHELFICEconserve,PARM_BOOL);
   
      
    %%% This code constructs the shelf ice load anomaly. Essentially the
    %%% model needs to know something about the pressure in the ocean at
    %%% the base of the ice shelf, and so we give it a "typical" pressure
    %%% field for that depth at each horizontal point. Once the model
    %%% starts running, it will continuously evolve that pressure, so the
    %%% pressure we prescribe here doesn't really matter - we just don't
    %%% want to blow up the model!
    salt_ref = 34;
    temp_ref = -1.9;   
    Pa1dbar = 1e4;
    SHELFICEloadAnomaly = zeros(Nx,Ny);
    for i=1:Nx      
      for j=1:Ny
        SHELFICEloadAnomaly(i,j) = 0;
        for k=1:Nr          
          if (zz(k) < icedraft(i,j))
            break;
          end
          Pressure = -rhoConst*g*zz(k);    
%           rhoShelfIce = densmdjwf(salt_ref,temp_ref,Pressure/Pa1dbar);
          rhoShelfIce = densmdjwf(sNorth(k),tNorth(k),Pressure/Pa1dbar);
          SHELFICEloadAnomaly(i,j) = SHELFICEloadAnomaly(i,j) + (g*(rhoShelfIce-rhoConst)*dz(k));                
         end
       end
    end
    
    figure(fignum);
    fignum = fignum+1;
    pcolor(X,Y,SHELFICEloadAnomaly);
    shading interp;
    colorbar;    
    
    %%% Write load anomaly 
    SHELFICEloadAnomalyFile='SHELFICEloadAnomalyFile.bin';
    writeDataset(SHELFICEloadAnomaly,fullfile(inputpath,SHELFICEloadAnomalyFile),ieee,prec);       
    shelfice_parm01.addParm('SHELFICEloadAnomalyFile','SHELFICEloadAnomalyFile.bin',PARM_STR); 
    
    %%% Ice shelf draft data
    writeDataset(icedraft,fullfile(inputpath,'SHELFICEtopoFile.bin'),ieee,prec);
    shelfice_parm01.addParm('SHELFICEtopoFile','SHELFICEtopoFile.bin',PARM_STR); 
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% WRITE THE 'data.shelfice' FILE %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    write_data_shelfice(inputpath,SHELFICE_PARM,listterm,realfmt);  

  end
  
  
  
  
  
    
    
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%% SEA ICE   %%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if (useSEAICE)
    % to store parameter names and values
    seaice_parm01 = parmlist;
    SEAICE_PARM = {seaice_parm01};

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%% SEA ICE  %%%%%%%%%%%%
      %%%%%%%% PARAMETERS %%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Oringinal albedos from llc_1/48th  other values from llc_2160
  % 

    SEAICEscaleSurfStress= true; % In the updated code (updated in Aug,2019),this issue has been solved. 20200121
                                 % By default, the sea ice stresses are not 
                                 % multiplied by the sea ice concentration. 
                                 % http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-August/011248.html
    SEAICEwriteState   = true;
    SEAICEuseDYNAMICS  = true;
    SEAICE_multDim     = 7;
    SEAICE_dryIceAlb   = 0.8783;
  %   SEAICE_dryIceAlb   = 0.8509;
    SEAICE_wetIceAlb   = 0.7869;
  %   SEAICE_wetIceAlb   = 0.7284;
    SEAICE_drySnowAlb  = 0.9482;
  %   SEAICE_drySnowAlb  = 0.7754;
    SEAICE_wetSnowAlb  = 0.8216;
  %   SEAICE_wetSnowAlb  = 0.7753;
    SEAICE_waterDrag   = 5.5399/1000; % water-ice drag coefficient (non-dim.)
    SEAICE_drag        = 0.002;   % air-ice drag coefficient (non-dim.)
    HO                 = 0.1; 
  %   HO                 = .05;

%     SEAICE_no_slip          = false;
    SEAICE_no_slip          = true;

  %   SEAICEadvScheme         = 7;
    SEAICEadvScheme         = 33;
    SEAICEmomAdvection      = false; % Default: false
    %%%SOSEdoesn't have a seaice dataset for salinity, they used this value
    %%%in their estimate

    LSR_ERROR               = 1.0e-5;  
    SEAICEnonLinIterMax     = 10;

    MIN_ATEMP               = -50;
    MIN_TICE                = -50;
    SEAICE_area_reg         = 0.15;
    SEAICE_hice_reg         = 0.1;
    IMAX_TICE               = 6;
    SEAICE_EPS		      = 1.0e-8;
  %   SEAICE_EPS              = 2.0e-9;
    SEAICE_doOpenWaterMelt  = true;
    SEAICE_areaLossFormula  = 1;
    SEAICE_wetAlbTemp       = 0.0;
    SEAICE_saltFrac         = 0.3;
  %   SEAICE_frazilFrac       = 0.003;
   SEAICE_frazilFrac       = 0.01;
  %   SEAICE_frazilFrac       = 1.0; % frazil to sea ice conversion rate, as fraction (relative to the local freezing point of sea ice water)


    % Initial fractional sea ice cover, range[0,1]; initializes variable AREA;
    Area = Ai0.*ones(Nx,Ny);
  %   Area(:,Ny-seaiceSpongeThickness:Ny) = 0;
    % Initial sea ice thickness averaged over grid cell in meters; initializes variable HEFF;
    Heff = Hi0.*ones(Nx,Ny); 
  %   Heff(:,Ny-seaiceSpongeThickness:Ny) = 0;
    % Initial snow thickness on sea ice averaged over grid cell in meters; initializes variable HSNOW;
    Hsnow = Hs0.*ones(Nx,Ny);
    % Initial salinity of sea ice averaged over grid cell in g/m^2; initializes variable HSALT;
    Hsalt = (Si0*rho_i*Hi0).*ones(Nx,Ny); 

    uIce = zeros(Nx,Ny); %%% Initial sea ice velosity
    vIce = zeros(Nx,Ny);

    AreaFile = 'AreaFile.bin';
    HeffFile = 'HeffFile.bin';
    HsnowFile = 'HsnowFile.bin';
    HsaltFile = 'HsaltFile.bin';
    uIceFile = 'uIceFile.bin';
    vIceFile = 'vIceFile.bin';  

    writeDataset(Area,fullfile(inputpath,AreaFile),ieee,prec);
    writeDataset(Heff,fullfile(inputpath,HeffFile),ieee,prec);
    writeDataset(Hsnow,fullfile(inputpath,HsnowFile),ieee,prec);
    writeDataset(Hsalt,fullfile(inputpath,HsaltFile),ieee,prec);  
    writeDataset(uIce,fullfile(inputpath,uIceFile),ieee,prec);
    writeDataset(vIce,fullfile(inputpath,vIceFile),ieee,prec); 

    seaice_parm01.addParm('SEAICEscaleSurfStress',SEAICEscaleSurfStress,PARM_BOOL);
    seaice_parm01.addParm('LSR_ERROR',LSR_ERROR,PARM_REAL);
    seaice_parm01.addParm('SEAICEnonLinIterMax',SEAICEnonLinIterMax,PARM_INT);
    seaice_parm01.addParm('SEAICEwriteState',SEAICEwriteState,PARM_BOOL);
    seaice_parm01.addParm('SEAICEuseDYNAMICS',SEAICEuseDYNAMICS,PARM_BOOL);
    seaice_parm01.addParm('SEAICE_multDim',SEAICE_multDim,PARM_INT);
    seaice_parm01.addParm('SEAICE_dryIceAlb',SEAICE_dryIceAlb,PARM_REAL);
    seaice_parm01.addParm('SEAICE_wetIceAlb',SEAICE_wetIceAlb,PARM_REAL);
    seaice_parm01.addParm('SEAICE_drySnowAlb',SEAICE_drySnowAlb,PARM_REAL);
    seaice_parm01.addParm('SEAICE_wetSnowAlb',SEAICE_wetSnowAlb,PARM_REAL);
    seaice_parm01.addParm('SEAICE_waterDrag',SEAICE_waterDrag,PARM_REAL);
    seaice_parm01.addParm('SEAICE_drag',SEAICE_drag,PARM_REAL);
    seaice_parm01.addParm('HO',HO,PARM_REAL);
  %   seaice_parm01.addParm('SEAICE_dryIceAlb_south',SEAICE_dryIceAlb_south,PARM_REAL);
  %   seaice_parm01.addParm('SEAICE_wetIceAlb_south',SEAICE_wetIceAlb_south,PARM_REAL);
  %   seaice_parm01.addParm('SEAICE_drySnowAlb_south',SEAICE_drySnowAlb_south,PARM_REAL);
  %   seaice_parm01.addParm('SEAICE_wetSnowAlb_south',SEAICE_wetSnowAlb_south,PARM_REAL);
  %   seaice_parm01.addParm('SEAICE_waterDrag_south',SEAICE_waterDrag_south,PARM_REAL);
  %   seaice_parm01.addParm('SEAICE_drag_south',SEAICE_drag_south,PARM_REAL);
    seaice_parm01.addParm('SEAICE_no_slip',SEAICE_no_slip,PARM_BOOL);
  %   seaice_parm01.addParm('SEAICE_salinity',SEAICE_salinity,PARM_REAL);
    seaice_parm01.addParm('SEAICEadvScheme',SEAICEadvScheme,PARM_INT);
    seaice_parm01.addParm('SEAICEmomAdvection',SEAICEmomAdvection,PARM_BOOL);
    seaice_parm01.addParm('MIN_ATEMP',MIN_ATEMP,PARM_REAL);
    seaice_parm01.addParm('MIN_TICE',MIN_TICE,PARM_REAL);
    seaice_parm01.addParm('SEAICE_area_reg',SEAICE_area_reg,PARM_REAL);
    seaice_parm01.addParm('SEAICE_hice_reg',SEAICE_hice_reg,PARM_REAL);
    seaice_parm01.addParm('IMAX_TICE',IMAX_TICE,PARM_INT);
    seaice_parm01.addParm('SEAICE_EPS',SEAICE_EPS,PARM_REAL);
    seaice_parm01.addParm('SEAICE_doOpenWaterMelt',SEAICE_doOpenWaterMelt,PARM_BOOL);
    seaice_parm01.addParm('SEAICE_areaLossFormula',SEAICE_areaLossFormula,PARM_INT);
    seaice_parm01.addParm('SEAICE_wetAlbTemp',SEAICE_wetAlbTemp,PARM_REAL);
    seaice_parm01.addParm('SEAICE_saltFrac',SEAICE_saltFrac,PARM_REAL);
    seaice_parm01.addParm('SEAICE_frazilFrac',SEAICE_frazilFrac,PARM_REAL);

    seaice_parm01.addParm('HeffFile',HeffFile,PARM_STR);
    seaice_parm01.addParm('AreaFile',AreaFile,PARM_STR);
    seaice_parm01.addParm('HsnowFile',HsnowFile,PARM_STR);
    seaice_parm01.addParm('HsaltFile',HsaltFile,PARM_STR);
    seaice_parm01.addParm('uIceFile',uIceFile,PARM_STR);
    seaice_parm01.addParm('vIceFile',vIceFile,PARM_STR);

  %   SEAICE_dalton = 0; % ice-ocean transfer coefficient for latent and sensible heat (non-dim.)
  %   seaice_parm01.addParm('SEAICE_dalton',SEAICE_dalton,PARM_REAL);
  %   SEAICErestoreUnderIce = true;  % enable restoring to climatology under ice
  %   seaice_parm01.addParm('SEAICErestoreUnderIce',SEAICErestoreUnderIce,PARM_BOOL);




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% WRITE THE 'data.seaice' FILE %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    write_data_seaice(inputpath,SEAICE_PARM,listterm,realfmt);  
    
  end
  
  
  
  
  



  %%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%
  %%%%%%EXF PKG%%%%%
  %%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%

  if (useEXF)
    
    %%% To store parameter names and values
    EXF_NML_01 = parmlist;
    EXF_NML_02 = parmlist;
    EXF_NML_03 = parmlist;
    EXF_NML_04 = parmlist;
    EXF_NML_OBCS = parmlist;
    EXF_PARM = {EXF_NML_01,EXF_NML_02,EXF_NML_03,EXF_NML_04,EXF_NML_OBCS}; 



    if(useSEAICE)

        exf_albedo = 0.15; 
        exf_scal_BulkCdn  = 1.015;
        exf_iprec         = 64;  
        useExfYearlyFields= false;
        useExfCheckRange  = false;
      %  	useRelativeWind   = true;
        useRelativeWind   = false;
          repeatPeriod      = 20*t1year;

      %     exf_offset_atemp =  273.16;
          %%%runoff from ERA is in hours, need to convert to seconds
      %     exf_inscal_runoff = 1.14e-04;

        EXF_NML_01.addParm('exf_albedo',exf_albedo,PARM_INT);
        if (EXFoption ~= 5)
           EXF_NML_01.addParm('exf_scal_BulkCdn',exf_scal_BulkCdn,PARM_REAL);
        end
        EXF_NML_01.addParm('exf_iprec',exf_iprec,PARM_INT);
        EXF_NML_01.addParm('useExfYearlyFields',useExfYearlyFields,PARM_BOOL);
        EXF_NML_01.addParm('useExfCheckRange',useExfCheckRange,PARM_BOOL);
        if(~useEXFwindstress)
            EXF_NML_01.addParm('useRelativeWind',useRelativeWind,PARM_BOOL);
        end
        EXF_NML_01.addParm('repeatPeriod',repeatPeriod,PARM_REAL);
      %   EXF_NML_03.addParm('exf_offset_atemp',exf_offset_atemp,PARM_REAL);
      %   EXF_NML_03.addParm('exf_inscal_runoff',exf_inscal_runoff,PARM_REAL);
        if (useEXFwindstress)
           readStressOnCgrid = true;
           EXF_NML_01.addParm('readStressOnCgrid',readStressOnCgrid,PARM_BOOL);
        end

          rho_a = 1.3;               %%% Air density, kg/m^3
      %     Ua = -6;
      %     Va = 6;

        if (Ua~=0)
      %     uwind = -sqrt(abs(tau_zonal)/rho_a/SEAICE_drag).*ones(Nx,Ny); % Zonal 10-m wind speed 
          uwind = [Ua:-Ua/(Ny-1):0].*ones(Nx,1); 
        else
          uwind = zeros(Nx,Ny); 
        end
        if (Va~=0)
      %     vwind = sqrt(abs(tau_merid)/rho_a/SEAICE_drag).*ones(Nx,Ny); % Meridional 10-m wind speed
           vwind = [Va:-Va/(Ny-1):0].*ones(Nx,1); 
        else
          vwind = zeros(Nx,Ny); 
        end

         %%% Plot the wind speed 
        if (showplots)
          figure(fignum);
          fignum = fignum + 1;
          clf;
          plot(yy/1000,uwind(1,:),'LineWidth',1.5);
          xlabel('Offshore distance (km)');
          ylabel('u_a');
          title('Zonal wind velocity (m/s)');
          set(gca,'fontsize',fontsize-1);
          PLOT = gcf;
          PLOT.Position = [263 149 567 336];  
        end    

         %%% Plot the wind speed 
        if (showplots)
          figure(fignum);
          fignum = fignum + 1;
          clf;
          plot(yy/1000,vwind(1,:),'LineWidth',1.5);
          xlabel('Offshore distance (km)');
          ylabel('v_a');
          title('Meridional wind velocity (m/s)');
          set(gca,'fontsize',fontsize-1);
          PLOT = gcf;
          PLOT.Position = [263 149 567 336];  
        end    


      %     Ur = sqrt((abs(tau_zonal)+abs(tau_merid))./rho_a./SEAICE_drag).*ones(Nx,Ny);
        uwindfile = 'uwindfile.bin';
        vwindfile = 'vwindfile.bin';
        writeDataset(uwind,fullfile(inputpath,uwindfile),ieee,prec);
        writeDataset(vwind,fullfile(inputpath,vwindfile),ieee,prec);
        EXF_NML_02.addParm('uwindfile',uwindfile,PARM_STR);
        EXF_NML_02.addParm('vwindfile',vwindfile,PARM_STR);  


      % Read-in atemp, aqh, swdown, lwdown, precip, and runoff. Compute hflux, swflux and sflux.
          Kice = 2.1656; %%% Ice thermal conductivity, W/(m*degK)
          ice_abs = 1-SEAICE_dryIceAlb; %%% Ice absorption

      % meanLWdown = zeros(1,size([-50:10],2));
      %  for  TaDegC = -50:10
          TaDegC = -10;
          Ta = 273.16+TaDegC; %%% Surface air temperature, degK
          Tw = 273.16+double(tNorth(1)); %%% surface water temperature
          TisDegC = -0.65;
      %     TisDegC = double(tNorth(1)); 
          Tis = 273.16+TisDegC; %%% Ice surface temperature
          SEAICE_emissivity = 0.970018; %%% Ice emissivity
          ocn_e = 5.50e-8 / 5.670e-8;  %%% Ocean emissivity, 0.97
          sigma = 5.67/10^8; %%% Stefan-Boltzman'n constant
          exf_iceCh = 1.63e-3; %%% sensible heat transfer coeff. over sea-ice   
          Cp_air = 1004; %%% Heat capacity at constant pressure 1004 J K-1 kg-1
          atemp = Ta.*ones(Nx,Ny); % Surface (2-m) air temperature in deg K
          aqh = 6.1094/(rho_a*287*Tis/100)*exp(17.625*TisDegC/(TisDegC+243.04)).*ones(Nx,Ny); % 0.0057, Surface (2m) specific humidity in kg/kg. Typical range: 0 < aqh < 0.02




      if(EXFoption == 3)
          swdown = 0.*ones(Nx,Ny); 
          precip = 0.*ones(Nx,Ny); 
          runoff = 0.*ones(Nx,Ny);  
      %     if(usezonalwind || useoffshorewind)
      %         SH = rho_a*Cp_air*exf_iceCh.*Ur.*abs(Ta-Tis);
      %         LH = -40.*SH/mean(mean(SH));
      %         lwdown = (ice_e*sigma*Tis^4 + SH + LH)/ice_e;
      %         SH = mean(mean(SH)).*ones(Nx,Ny);
      %         LH = mean(mean(LH)).*ones(Nx,Ny);
      %         lwdown =  mean(mean(lwdown)).*ones(Nx,Ny);
      %         meanSHdown = mean(mean(SH))
      %         meanLHdown = mean(mean(LH))
      %     else
      %         lwdown = (ice_e*sigma*Tis^4)/ice_e.*ones(Nx,Ny);
      %     end
      %     meanLWdown = mean(mean(lwdown))
      %     lwdown = 320.*ones(Nx,Ny);  

          CondHeat = 0.5*1/Hi0;  %%% SItice ~ -1.62 degC, Tio ~ -1.87 degC => Conductive heat flux from ice surface to ocean is about 0.5 W/m^2
          lwdown = (CondHeat/ice_abs+320)*ones(Nx,Ny);
          lwdown(1)

          if (Hi0==0)
              lwdown = 324.1085*ones(Nx,Ny);
              lwdown(1)
          end



        %%% Plot the downward longwave radiation in W/m^2
        if (showplots)
          figure(fignum);
          fignum = fignum + 1;
          clf;
          plot(yy/1000,squeeze(lwdown(1,:)),'LineWidth',1.5);
          xlabel('Offshore distance (km)');
          ylabel('LWdown (W/m^2)');
          title('Downward longwave radiation');
          set(gca,'fontsize',fontsize-1);
          PLOT = gcf;
          PLOT.Position = [263 149 567 336];  
        end
          %%% Save the figure
          savefig([imgpath '/LWdown.fig']);
          saveas(gcf,[imgpath '/LWdown.png']);

          atempfile  = 'atempfile.bin';
          aqhfile    = 'aqhfile.bin';
          swdownfile = 'swdownfile.bin';
          lwdownfile = 'lwdownfile.bin';
          precipfile = 'precipfile.bin'; 
          runofffile = 'runofffile.bin';
          writeDataset(atemp,fullfile(inputpath,atempfile),ieee,prec);
          writeDataset(aqh,fullfile(inputpath,aqhfile),ieee,prec);
          writeDataset(swdown,fullfile(inputpath,swdownfile),ieee,prec);
          writeDataset(lwdown,fullfile(inputpath,lwdownfile),ieee,prec);
          writeDataset(precip,fullfile(inputpath,precipfile),ieee,prec);
          writeDataset(runoff,fullfile(inputpath,runofffile),ieee,prec);
          EXF_NML_02.addParm('atempfile',atempfile,PARM_STR);
          EXF_NML_02.addParm('aqhfile',aqhfile,PARM_STR);
          EXF_NML_02.addParm('swdownfile',swdownfile,PARM_STR);
          EXF_NML_02.addParm('lwdownfile',lwdownfile,PARM_STR);
          EXF_NML_02.addParm('precipfile',precipfile,PARM_STR);
          EXF_NML_02.addParm('runofffile',runofffile,PARM_STR);

      elseif(EXFoption == 5)    
          hs = 0.*ones(Nx,Ny); 
          hl = 0.*ones(Nx,Ny);  
          swdown = 0.*ones(Nx,Ny); 
          Tisurf = 273.16+temp_relax(:,:,1); %%% assume ice temperature at surface equals bottom
      %     lwdown = SEAICE_emissivity*sigma.*Tisurf.^4/(1-ice_abs);
          lwdown = SEAICE_emissivity*sigma.*Tisurf.^4;
          evap = 0.*ones(Nx,Ny);  
          precip = 0.*ones(Nx,Ny); 
          runoff = 0.*ones(Nx,Ny);  
          hsfile  = 'hsfile.bin';
          hlfile  = 'hlfile.bin';
          swdownfile = 'swdownfile.bin';
          lwdownfile = 'lwdownfile.bin';
          evapfile  = 'evapfile.bin';
          precipfile = 'precipfile.bin'; 
          runofffile = 'runofffile.bin'; 
          writeDataset(hs,fullfile(inputpath,hsfile),ieee,prec);
          writeDataset(hl,fullfile(inputpath,hlfile),ieee,prec);
          writeDataset(swdown,fullfile(inputpath,swdownfile),ieee,prec);
          writeDataset(lwdown,fullfile(inputpath,lwdownfile),ieee,prec);
          writeDataset(evap,fullfile(inputpath,evapfile),ieee,prec);    
          writeDataset(precip,fullfile(inputpath,precipfile),ieee,prec);
          writeDataset(runoff,fullfile(inputpath,runofffile),ieee,prec);
          EXF_NML_02.addParm('hsfile',hsfile,PARM_STR);
          EXF_NML_02.addParm('hlfile',hlfile,PARM_STR);
          EXF_NML_02.addParm('swdownfile',swdownfile,PARM_STR);
          EXF_NML_02.addParm('lwdownfile',lwdownfile,PARM_STR);
          EXF_NML_02.addParm('evapfile',evapfile,PARM_STR); 
          EXF_NML_02.addParm('precipfile',precipfile,PARM_STR);
          EXF_NML_02.addParm('runofffile',runofffile,PARM_STR);


      elseif (EXFoption == 6)
          hs = 0.*ones(Nx,Ny); 
          hl = 0.*ones(Nx,Ny);  
          swflux = 0.*ones(Nx,Ny);  
          lwflux = 0.*ones(Nx,Ny);  
          evap = 0.*ones(Nx,Ny);  
          precip = 0.*ones(Nx,Ny); 
          runoff = 0.*ones(Nx,Ny);  
          hsfile  = 'hsfile.bin';
          hlfile  = 'hlfile.bin';
          swfluxfile  = 'swfluxfile.bin';
          lwfluxfile  = 'lwfluxfile.bin';
          evapfile  = 'evapfile.bin';
          precipfile = 'precipfile.bin'; 
          runofffile = 'runofffile.bin'; 
          writeDataset(hs,fullfile(inputpath,hsfile),ieee,prec);
          writeDataset(hl,fullfile(inputpath,hlfile),ieee,prec);
          writeDataset(swflux,fullfile(inputpath,swfluxfile),ieee,prec);
          writeDataset(lwflux,fullfile(inputpath,lwfluxfile),ieee,prec);
          writeDataset(evap,fullfile(inputpath,evapfile),ieee,prec);    
          writeDataset(precip,fullfile(inputpath,precipfile),ieee,prec);
          writeDataset(runoff,fullfile(inputpath,runofffile),ieee,prec);
          EXF_NML_02.addParm('hsfile',hsfile,PARM_STR);
          EXF_NML_02.addParm('hlfile',hlfile,PARM_STR);
          EXF_NML_02.addParm('swfluxfile',swfluxfile,PARM_STR);
          EXF_NML_02.addParm('lwfluxfile',lwfluxfile,PARM_STR);
          EXF_NML_02.addParm('evapfile',evapfile,PARM_STR); 
          EXF_NML_02.addParm('precipfile',precipfile,PARM_STR);
          EXF_NML_02.addParm('runofffile',runofffile,PARM_STR);

      elseif (EXFoption == 4)   
          swflux = 0.*ones(Nx,Ny);  
          lwflux = 0.*ones(Nx,Ny);  
          precip = 0.*ones(Nx,Ny); 
          runoff = 0.*ones(Nx,Ny);  
          atempfile  = 'atempfile.bin';
          aqhfile    = 'aqhfile.bin';
          swfluxfile  = 'swfluxfile.bin';
          lwfluxfile  = 'lwfluxfile.bin';
          precipfile = 'precipfile.bin'; 
          runofffile = 'runofffile.bin'; 
          writeDataset(atemp,fullfile(inputpath,atempfile),ieee,prec);
          writeDataset(aqh,fullfile(inputpath,aqhfile),ieee,prec);
          writeDataset(swflux,fullfile(inputpath,swfluxfile),ieee,prec);
          writeDataset(lwflux,fullfile(inputpath,lwfluxfile),ieee,prec);
          writeDataset(precip,fullfile(inputpath,precipfile),ieee,prec);
          writeDataset(runoff,fullfile(inputpath,runofffile),ieee,prec);
          EXF_NML_02.addParm('atempfile',atempfile,PARM_STR);
          EXF_NML_02.addParm('aqhfile',aqhfile,PARM_STR);
          EXF_NML_02.addParm('swfluxfile',swfluxfile,PARM_STR);
          EXF_NML_02.addParm('lwfluxfile',lwfluxfile,PARM_STR);
          EXF_NML_02.addParm('precipfile',precipfile,PARM_STR);
          EXF_NML_02.addParm('runofffile',runofffile,PARM_STR);

      end

    end


    %%% Create the data.exf file
    write_data_exf(inputpath,EXF_PARM,listterm,realfmt);

  end
  
  
 
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% LAYERS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  if (useLAYERS)

    %%% To store parameter names and values
    layers_parm01 = parmlist;
    LAYERS_PARM = {layers_parm01};


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% LAYERS PARAMETERS %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%% Define parameters for layers package %%%

    %%% Number of fields for which to calculate layer fluxes  
    layers_maxNum = 1;

    %%% Specify potential density
    layers_name = char('RHO'); 

    layers_bounds = [0 35 ...
           35.8:0.05:36.3 ...
           36.4 36.54:0.02:36.66 ...
           36.7 36.73 36.76 36.8:0.1:37.1 ...
           37.13:0.02:37.17 37.18:0.004:37.206 ...
           37.21:0.003:37.3 37.5 40]; 


    %%% Reference level for calculation of potential density  
    layers_krho = 51; % High-resolution zz (51)=-1.9943e+03 m;

    %%% If set true, the GM bolus velocity is added to the calculation
    layers_bolus = false;  

    %%% Layers
    layers_parm01.addParm(['layers_bounds'],layers_bounds,PARM_REALS); 
    layers_parm01.addParm(['layers_krho'],layers_krho,PARM_INT); 
    layers_parm01.addParm(['layers_name'],strtrim(layers_name),PARM_STR); 
    layers_parm01.addParm('layers_bolus',layers_bolus,PARM_BOOL); 

    %%z% Create the data.layers file
    write_data_layers(inputpath,LAYERS_PARM,listterm,realfmt);

    %%% Create the LAYERS_SIZE.h file
    createLAYERSSIZEh(codepath,length(layers_bounds)-1,layers_maxNum); 

  end
  
  

 
 
 
 
 
 
 
 
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   
  %%% To store parameter names and values
  diag_parm01 = parmlist;
  diag_parm02 = parmlist;
  DIAG_PARM = {diag_parm01,diag_parm02};
  diag_matlab_parm01 = parmlist;
  DIAG_MATLAB_PARM = {diag_matlab_parm01}; %%% Matlab parameters need to be different
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  %%% Stores total number of diagnostic quantities
  ndiags = 0;
       
  diag_parm01.addParm('diag_mnc',false,PARM_BOOL);  
  diag_parm01.addParm('diag_pickup_read',false,PARM_BOOL);  
  diag_parm01.addParm('diag_pickup_write',false,PARM_BOOL);  
%   diag_parm01.addParm('diag_pickup_read_mnc',false,PARM_BOOL);  
%   diag_parm01.addParm('diag_pickup_write_mnc',false,PARM_BOOL); 
   
        


%          'LaTh1RHO','LaTz1RHO','LTha1RHO','LTza1RHO','LTto1RHO','LaSh1RHO','LaSz1RHO','LSha1RHO','LSza1RHO','LSto1RHO',...
%          'LaTh2TH','LaTz2TH','LTha2TH','LTza2TH','LTto2TH'...


%%% Annual mean diagnostics
diag_fields_avg = {...   
% % %      'UVEL','VVEL','WVEL','SALT','THETA','PHIHYD','ETAN',...%%% Basic state 
% % %      'TOTTTEND','TFLUX','ADVy_TH','VVELTH','oceQnet',...%%% Heat budget
% % %      'UVELSQ','VVELSQ','WVELSQ','UV_VEL_Z','WU_VEL','WV_VEL',...%%% Energy budget
% % %      'LaVH1RHO','LaHs1RHO',...%%% Overturning circ
% % %      'oceTAUX','oceTAUY','Um_dPhiX','Um_Advec','Um_Diss','Um_Ext',...
% % %      'Um_Cori','Um_AdvZ3','Um_AdvRe',...%%% Momentum budget
% % %      'SIarea','SIheff','SIuice','SItices','SIvice','SItaux','SItauy','SIsig12',...
% % %      'SIempmr','oceSflux','SIatmTx','SIatmTy','SIqnet','SIatmQnt',...%%% Sea ice
% % %       ...
%     ... %%%%%%%%% for spin-up
%     'TOTTTEND','TFLUX','VVELTH','ADVy_TH','oceQnet',...
%     'SIarea','SIheff','SIuice','SIvice','SIempmr','oceSflux',...
%     'UVELSQ','VVELSQ','WVELSQ'...
% % % % % %       ... %%%%%%%%% for analysis
% % % % % %       ... %%% Heat budget
% % % % % %          'TOTTTEND','TFLUX','KPPg_TH','oceQsw','WTHMASS',...
% % % % % %          'ADVr_TH','ADVx_TH','ADVy_TH','DFxE_TH','DFyE_TH','DFrI_TH','DFrE_TH',...
% % % % % %          ...
% % % % % %          'VVELTH', ...
% % % % % %          'oceQnet','UVELTH','WVELTH',...
% % % % % %       ... %%% Energy budget
% % % % % %          'UVELSQ','VVELSQ','WVELSQ',...
% % % % % %          'UV_VEL_Z','WU_VEL','WV_VEL',...
% % % % % %       ... %%% Sea ice
% % % % % %          'SIarea','SIheff','SIuice','SIvice','SIsig12',...
% % % % % %          'SItices','SIqnet','SIempmr','SIatmQnt',...
% % % % % %          ...
% % % % % %       ... %%% Salt budget
% % % % % %          'TOTSTEND','SFLUX','KPPg_SLT','oceFWflx','WSLTMASS',...
% % % % % %          'ADVr_SLT','ADVx_SLT','ADVy_SLT','DFrE_SLT','DFxE_SLT','DFyE_SLT','DFrI_SLT',...
% % % % % %          ...
% % % % % %          'VVELSLT',...
% % % % % %          'oceSflux','UVELSLT','WVELSLT',...
% % % % % %       ... %%% Momentum budget
% % % % % %          'ETAN',...
% % % % % %          'oceTAUX','oceTAUY',...
% % % % % %      ... %%% Overturning streamfunction
% % % % % %          'RHOAnoma','LaUH1RHO','LaHw1RHO','LaTr1RHO',... 
% % % % % %                     'LaUH2TH','LaHw2TH',... 
% % % % % %      ...
% % % % % %          'Um_Diss','Um_Advec','Um_dPhiX','Um_Ext',...
% % % % % %          'SItaux','SItauy','SIatmTx','SIatmTy',...   
% % % % % %          'Vm_Diss','Vm_Advec','Vm_Cori','Vm_dPhiY','Vm_Ext','Vm_AdvZ3','Vm_AdvRe',...
% % % % % %          'VISrI_Um','VISrI_Vm',...
     };
      
  numdiags_avg = length(diag_fields_avg);  
  diag_freq_avg = 1*t1year;

  diag_phase_avg = 0;    
      
  for n=1:numdiags_avg    
    ndiags = ndiags + 1;
    diag_parm01.addParm(['fields(1,',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields{1,',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames{',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL);   
  end
  
  
  %%%%%% Daily output 
  diag_fields_inst = {...
    'UVEL','VVEL','WVEL','THETA','SALT','ETAN', ...
    %'SIarea','SIheff','SIuice','SIvice' ...
      };
  numdiags_inst = length(diag_fields_inst);  
  diag_freq_inst =1*t1year/12;
  diag_phase_inst = 0;
  
  for n=1:numdiags_inst    
    ndiags = ndiags + 1;
    diag_parm01.addParm(['fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL);     
  end
  



  %%% Create the data.diagnostics file
  write_data_diagnostics(inputpath,DIAG_PARM,listterm,realfmt);
  
  %%% Create the DIAGNOSTICS_SIZE.h file
  if(useLAYERS)
    createDIAGSIZEh(codepath,ndiags,max(Nr,length(layers_bounds)-1));
  else
    createDIAGSIZEh(codepath,ndiags,Nr);
  end
  
  
  
  
  



  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  %%%%% OBCS %%%%%
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  
  %%% TODO configure with a single OB at the northern end of the domain
  
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% OBCS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  %%% Add 2-element cell arrays to this cell array in the form 
  %%%  OBCS_PARM{1} = addParameter(OBCS_PARM{1},'paramName',paramValue,parmType);
  %%% to specify additional parameters. The parameter type parmType must
  %%% take one of the integer values above.    
  obcs_parm01 = parmlist;
  obcs_parm02 = parmlist;
  obcs_parm03 = parmlist;
  obcs_parm04 = parmlist;
  obcs_parm05 = parmlist;
  OBCS_PARM = {obcs_parm01,obcs_parm02,obcs_parm03,obcs_parm04,obcs_parm05};  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DEFINE OPEN BOUNDARY TYPES (OBCS_PARM01) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%% Set boundary points that are open     
  OB_Jnorth= -1*ones(1,Nx);
  obcs_parm01.addParm('OB_Jnorth',OB_Jnorth,PARM_INTS);    
  
  %%% Enables/disables Orlanski radiation conditions at the boundaries -
  %%% allows waves to propagate out through the boundary with minimal
  %%% reflection      
  obcs_parm01.addParm('useOrlanskiNorth',useOrlanskiNorth,PARM_BOOL);  
  
  %%% Enforces mass conservation across the northern boundary by adding a
  %%% barotropic inflow/outflow    
  obcs_parm01.addParm('useOBCSbalance',useOBCSbalance,PARM_BOOL);
  if (useOBCSbalance)      
      OBCS_balanceFacN = -1; %%% A value -1 balances an individual boundary       
      obcs_parm01.addParm('OBCS_balanceFacN',OBCS_balanceFacN,PARM_REAL);         
  end
  
  %%% Enables/disables sponge layers   
  useOBCSsponge = true;
  obcs_parm01.addParm('useOBCSsponge',useOBCSsponge,PARM_BOOL);
    
  if(useSEAICE)
      useSeaiceSponge = true;
      obcs_parm01.addParm('useSeaiceSponge',useSeaiceSponge,PARM_BOOL);
  else 
      useSeaiceSponge = false;
  end
  
  %%% Set boundary velocities and temperatures
  useOBCSprescribe = true;  
  
  %%% Set northern boundary properties
  OBNt = ones(Nx,1)*tNorth;
  OBNs = ones(Nx,1)*sNorth;  
  
  %%% Write boundary variables to files  
  writeDataset(OBNt,fullfile(inputpath,'OBNtFile.bin'),ieee,prec);
  writeDataset(OBNs,fullfile(inputpath,'OBNsFile.bin'),ieee,prec);  

  %%% Set OBCS prescription parameters
  obcs_parm01.addParm('useOBCSprescribe',useOBCSprescribe,PARM_BOOL);
  obcs_parm01.addParm('OBNtFile','OBNtFile.bin',PARM_STR);
  obcs_parm01.addParm('OBNsFile','OBNsFile.bin',PARM_STR);  

  %%% TODO
    if(useSEAICE)
%     OBNa = Ai0.*ones(Nx,1);
%     OBNh = Hi0.*ones(Nx,1);
%     OBNsn = Hs0.*ones(Nx,1); %%% snow thickness
%     OBNsl = Si0.*ones(Nx,1); %%% sea ice salinity
%     OBNuice = -0.24*ones(Nx,1);
%     OBNuice = 0*ones(Nx,1);
%     OBNvice = 0.005.*ones(Nx,1);
%     writeDataset(OBNa,fullfile(inputpath,'OBNaFile.bin'),ieee,prec);
%     writeDataset(OBNh,fullfile(inputpath,'OBNhFile.bin'),ieee,prec);
%     writeDataset(OBNsn,fullfile(inputpath,'OBNsnFile.bin'),ieee,prec);
%     writeDataset(OBNsl,fullfile(inputpath,'OBNslFile.bin'),ieee,prec);
%     writeDataset(OBNuice,fullfile(inputpath,'OBNuiceFile.bin'),ieee,prec);
%     writeDataset(OBNvice,fullfile(inputpath,'OBNviceFile.bin'),ieee,prec);
%     obcs_parm01.addParm('OBNaFile','OBNaFile.bin',PARM_STR);
%     obcs_parm01.addParm('OBNhFile','OBNhFile.bin',PARM_STR);
%     obcs_parm01.addParm('OBNsnFile','OBNsnFile.bin',PARM_STR);
%     obcs_parm01.addParm('OBNslFile','OBNslFile.bin',PARM_STR);
%     obcs_parm01.addParm('OBNuiceFile','OBNuiceFile.bin',PARM_STR);
%     obcs_parm01.addParm('OBNviceFile','OBNviceFile.bin',PARM_STR);

%%% Calculate free-drift ice velocities at the southern boundary, ignoring
%%% ice internal stress, pressure caused by sea surface hight variation.
%%% Assume zonal ocean velocity at the coast = 0, 
%%% meridional ocean at the coast averaged over tidal cycles = 0,            
%%% fractional ice cover Ai0 = 1.
    rho_o = 1027;         %%% Water density, kg/m^3
    tao_aix = rho_a*SEAICE_drag*sqrt(Ua^2+Va^2)*Ua;       %%% Air-ice stress in x direction, N/m2
    tao_aiy = rho_a*SEAICE_drag*sqrt(Ua^2+Va^2)*Va;       %%% Air-ice stress in y direction, N/m2
    syms ui vi
    eq1 =  rho_i*Hi0*f0*vi + tao_aix - rho_o*SEAICE_waterDrag*sqrt(ui^2+vi^2)*ui;
    eq2 = -rho_i*Hi0*f0*ui + tao_aiy - rho_o*SEAICE_waterDrag*sqrt(ui^2+vi^2)*vi;
    eqns = [eq1, eq2];
    [solui solvi] = solve(eqns,[ui vi]);
    Sui = double(real(solui));
    Svi = double(real(solvi));
    ui_idx = (Sui<0);
    obsuice = Sui(ui_idx)
    obsvice = Svi(ui_idx)
    

    if (useRelaxSurfT)
        obsuice = 0
        obsvice = 0
    end
   
    if(Ua == 0 && Va == 0)
        obsuice = 0
        obsvice = 0
    end
    
    OBSa = Ai0.*ones(Nx,1);
    OBSh = Hi0.*ones(Nx,1);
    OBSsn = Hs0.*ones(Nx,1); %%% snow thickness
    OBSsl = Si0.*ones(Nx,1); %%% sea ice salinity
    OBSuice = obsuice.*ones(Nx,1); %%% Initial zonal ice velocity should be westward (negative!) or zero.
    OBSvice = obsvice.*ones(Nx,1);
    writeDataset(OBSa,fullfile(inputpath,'OBSaFile.bin'),ieee,prec);
    writeDataset(OBSh,fullfile(inputpath,'OBShFile.bin'),ieee,prec);
    writeDataset(OBSsn,fullfile(inputpath,'OBSsnFile.bin'),ieee,prec);
    writeDataset(OBSsl,fullfile(inputpath,'OBSslFile.bin'),ieee,prec);
    writeDataset(OBSuice,fullfile(inputpath,'OBSuiceFile.bin'),ieee,prec);
    writeDataset(OBSvice,fullfile(inputpath,'OBSviceFile.bin'),ieee,prec);
    obcs_parm01.addParm('OBSaFile','OBSaFile.bin',PARM_STR);
    obcs_parm01.addParm('OBShFile','OBShFile.bin',PARM_STR);
    obcs_parm01.addParm('OBSsnFile','OBSsnFile.bin',PARM_STR);
    obcs_parm01.addParm('OBSslFile','OBSslFile.bin',PARM_STR);
    obcs_parm01.addParm('OBSuiceFile','OBSuiceFile.bin',PARM_STR);
    obcs_parm01.addParm('OBSviceFile','OBSviceFile.bin',PARM_STR);
    
    
    
    
    OBNa = Ai0.*ones(Nx,1);
    OBNh = Hi0.*ones(Nx,1);
    OBNsn = Hs0.*ones(Nx,1); %%% snow thickness
    OBNsl = Si0.*ones(Nx,1); %%% sea ice salinity
    OBNuice = obsuice.*ones(Nx,1); %%% Initial zonal ice velocity should be westward (negative!) or zero.
    OBNvice = obsvice.*ones(Nx,1);
    
    writeDataset(OBNa,fullfile(inputpath,'OBNaFile.bin'),ieee,prec);
    writeDataset(OBNh,fullfile(inputpath,'OBNhFile.bin'),ieee,prec);
    writeDataset(OBNsn,fullfile(inputpath,'OBNsnFile.bin'),ieee,prec);
    writeDataset(OBNsl,fullfile(inputpath,'OBNslFile.bin'),ieee,prec);
    writeDataset(OBNuice,fullfile(inputpath,'OBNuiceFile.bin'),ieee,prec);
    writeDataset(OBNvice,fullfile(inputpath,'OBNviceFile.bin'),ieee,prec);
    obcs_parm01.addParm('OBNaFile','OBNaFile.bin',PARM_STR);
    obcs_parm01.addParm('OBNhFile','OBNhFile.bin',PARM_STR);
    obcs_parm01.addParm('OBNsnFile','OBNsnFile.bin',PARM_STR);
    obcs_parm01.addParm('OBNslFile','OBNslFile.bin',PARM_STR);
    obcs_parm01.addParm('OBNuiceFile','OBNuiceFile.bin',PARM_STR);
    obcs_parm01.addParm('OBNviceFile','OBNviceFile.bin',PARM_STR);
    
    
    end

%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %%%%% ORLANSKI OPTIONS (OBCS_PARM02) %%%%%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   
%   %%% Velocity averaging time scale - must be larger than deltaT.
%   %%% The Orlanski radiation condition computes the characteristic velocity
%   %%% at the boundary by averaging the spatial derivative normal to the 
%   %%% boundary divided by the time step over this period.
%   %%% At the moment we're using the magic engineering factor of 3.
%   cvelTimeScale = 2*deltaT;
%   %%% Max dimensionless CFL for Adams-Bashforth 2nd-order method
%   CMAX = 0.45; 
%   
%   obcs_parm02.addParm('cvelTimeScale',cvelTimeScale,PARM_REAL);
%   obcs_parm02.addParm('CMAX',CMAX,PARM_REAL);
  
  


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% SPONGE LAYER OPTIONS (OBCS_PARM03) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  if (useOBCSsponge)

    obcs_parm03.addParm('spongeThickness',spongeThickness,PARM_INT);

    %%% Urelaxobcsbound/inner set the relaxation time scales for the Eastern and Western boundaries, 
    %%% Vrelaxobcsbound/inner for the Northern and Southern boundaries.
  
    %%% Relaxation time scales
    Vrelaxobcsinner = 864000;
    Vrelaxobcsbound = 43200; 
    obcs_parm03.addParm('Vrelaxobcsinner',Vrelaxobcsinner,PARM_REAL);
    obcs_parm03.addParm('Vrelaxobcsbound',Vrelaxobcsbound,PARM_REAL);
    
  end
    
  %%% TODO
  if (useSeaiceSponge)
       T_relaxinner = 864000/10;
       T_relaxbound = 43200/6;
    Arelaxobcsinner = T_relaxinner;
    Arelaxobcsbound = T_relaxbound;
    Hrelaxobcsinner = T_relaxinner;
    Hrelaxobcsbound = T_relaxbound;
    SLrelaxobcsinner = T_relaxinner;
    SLrelaxobcsbound = T_relaxbound;
    SNrelaxobcsinner = T_relaxinner;
    SNrelaxobcsbound = T_relaxbound;
    obcs_parm05.addParm('seaiceSpongeThickness',seaiceSpongeThickness,PARM_INT);
    obcs_parm05.addParm('Arelaxobcsinner',Arelaxobcsinner,PARM_REAL);
    obcs_parm05.addParm('Arelaxobcsbound',Arelaxobcsbound,PARM_REAL);
    obcs_parm05.addParm('Hrelaxobcsinner',Hrelaxobcsinner,PARM_REAL);
    obcs_parm05.addParm('Hrelaxobcsbound',Hrelaxobcsbound,PARM_REAL);
    obcs_parm05.addParm('SLrelaxobcsinner',SLrelaxobcsinner,PARM_REAL);
    obcs_parm05.addParm('SLrelaxobcsbound',SLrelaxobcsbound,PARM_REAL);
    obcs_parm05.addParm('SNrelaxobcsinner',SNrelaxobcsinner,PARM_REAL);
    obcs_parm05.addParm('SNrelaxobcsbound',SNrelaxobcsbound,PARM_REAL);
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.obcs' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  %%% Creates the 'data.obcs' file
  write_data_obcs(inputpath,OBCS_PARM,listterm,realfmt);

  
  
  
  
  
  
  
  

  
  %%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%
  %%%%% PACKAGES %%%%%
  %%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%
  
  packages = parmlist;
  PACKAGE_PARM = {packages};  
  
  packages.addParm('useDiagnostics',true,PARM_BOOL);    
  packages.addParm('useKPP',true,PARM_BOOL);  
  packages.addParm('useEXF',useEXF,PARM_BOOL);        
  packages.addParm('useCAL',useEXF,PARM_BOOL); 
  packages.addParm('useSHELFICE',useSHELFICE,PARM_BOOL);
  packages.addParm('useSEAICE',useSEAICE,PARM_BOOL);
  packages.addParm('useOBCS',useOBCS,PARM_BOOL);  
  packages.addParm('useLAYERS',useLAYERS,PARM_BOOL);  

  %%% Create the data.pkg file
  write_data_pkg(inputpath,PACKAGE_PARM,listterm,realfmt);
  
  
  
  
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE PARAMETERS TO A MATLAB FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%% Creates a matlab file defining all input parameters
  ALL_PARMS =[PARM PACKAGE_PARM DIAG_MATLAB_PARM];
  
  if (useEXF)
    ALL_PARMS = [ALL_PARMS EXF_PARM];
  end
  if (useSEAICE)
    ALL_PARMS = [ALL_PARMS SEAICE_PARM];
  end  
  if (useSHELFICE)
    ALL_PARMS = [ALL_PARMS SHELFICE_PARM];
  end  
  if (useLAYERS)
    ALL_PARMS = [ALL_PARMS LAYERS_PARM];
  end  
 
  %%% Creates a matlab file defining all input parameters
  write_matlab_params(inputpath,ALL_PARMS,realfmt);
  
  
end
  
