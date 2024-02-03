%%%
%%% anim3Dstrat.m
%%%
%%% Makes a 3D movie of our model stratification.
%%%
%%% NOTE: Doesn't account for u/v gridpoint locations, and doesn't handle
%%% partial cells.
%%%

%%% Select potential temperature surface
theta_plot = 0.4;

%%% Select simulation
%expdir = '/home/garrett/Projects/MITgcm_ISC/experiments/slope375-GLIB-explore-18';
%expname = 'at125';
%%% Read experiment data

loadexp;

%%% Vertical grid spacing matrix
DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);

%%% Diagnostic indix corresponding to instantaneous velocity
diagnum = length(diag_frequency);

%%% This needs to be set to ensure we are using the correct output
%%% frequency
diagfreq = diag_frequency(1);

%%% Frequency of diagnostic output
dumpFreq = abs(diagfreq);
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Plotting options
fontsize = 16;
figlabel = '';
M = moviein(nDumps);
  
%%% Read snapshot
for n=5:length(dumpIters)
% for n = 1

%   uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),dumpIters(n));      
%   vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),dumpIters(n));   
  theta = rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters(n));         
  tt =  (dumpIters(n)-dumpIters(1))*deltaT/86400

  %%% Load time-mean output
  %load(fullfile(prodir,[expname,'_tavg_5yrs.mat']),'PHIHYD');
  %pp = PHIHYD;
  %eta = pp(:,:,1)/gravity;
  %eta(:,1) = NaN;
  %eta(:,end) = NaN;
  %%eta = eta - mean(eta(:,end-1));

  %%% Remove topography
  theta(hFacC==0) = NaN;
  %eta(hFacC(:,:,1)==0) = NaN;

  fig = figure(1);
  clf;
  scrsz = get(0,'ScreenSize');
%   set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 600 600]);
  set(gcf,'Position',[780   129   786   666]);
  set(gcf,'Color','w');
  [YY,XX,ZZ]=meshgrid(yy,xx,zz);
  XX = XX / 1000;
  YY = YY / 1000;

%   %%% Calculate vorticity
%   ff = f0+beta*YY;
%   vort = zeros(Nx,Ny,Nr);
%   vort(:,1:Ny-1,:) = - (uvel(:,2:Ny,:)-uvel(:,1:Ny-1,:))/delY(1);
%   vort = vort + (vvel([2:Nx 1],:,:)-vvel(:,:,:))/delX(1);
%   vort(hFacS==0) = 0;
%   vort(hFacW==0) = 0;
%   vort(hFacC==0) = 0;
%   vort(:,Ny-1,:) = 0;
%   vort(:,2,:) = 0;
%   vort = vort ./ abs(ff);

  %%% Bathymetry
  [Y,X] = meshgrid(yy,xx);  
  % surf(X(:,2:end-1)/1000,Y(:,2:end-1)/1000,bathy(:,2:end-1));
  p = surface(X(:,2:end-1)/1000,Y(:,2:end-1)/1000,bathy(:,2:end-1));
  alpha(p,1)
  % p.FaceColor = [8*16+5 5*16+7 2*16+3]/255;
  p.FaceColor = [11*16+9 9*16+12 6*16+11]/255;
  p.EdgeColor = 'none';       

  icedraft=SHELFICEtopo
  icedraft_plot = icedraft;
  icedraft_plot(icedraft==0) = NaN;
  icetop_plot = 0*icedraft_plot;
  for i=1:Nx
    j = find(~isnan(icetop_plot(i,:)),1,'last');
    if (isempty(j))
      continue;
    else
      icetop_plot(i,j+1) = max(-200,bathy(i,j+1));
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


  hold on;

  %%% Plot SST
  %% p = surface(X(:,2:end-1)/1000,Y(:,2:end-1)/1000,0*X(:,2:end-1),theta(:,2:end-1,1));
  %% p.FaceColor = 'texturemap';
  %% colormap(pmkmp(56,'Swtth'));
  %% caxis([-2 12]);
  %% p.EdgeColor = 'none';         
  %% alpha(p,0.8);
  
%   %%% Plot Vorticity
%   p = surface(X(:,2:end-1)/1000,Y(:,2:end-1)/1000,0*X(:,2:end-1),vort(:,2:end-1,1));
%   p.FaceColor = 'texturemap';
%   colormap(redblue);
%   caxis([-2 12]);
%   p.EdgeColor = 'none';         
%   alpha(p,0.8);

  i_end = 1;
  theta_end = squeeze(theta(i_end,:,:));
  theta_end(theta_end==0) = NaN;
  [ZZZ,YYY] = meshgrid(zz,yy);
  zz_i = -H+5:10:-5;
  Nz_i = length(zz_i);
  [ZZZ_i,YYY_i] = meshgrid(zz_i,yy);
  theta_i = zeros(Ny,Nz_i);
  for j=1:Ny
    theta_i(j,:) = interp1(zz,theta_end(j,:),zz_i,'linear');
  end
  p = surface(xx(i_end)/1000*ones(size(YYY_i)),YYY_i/1000,ZZZ_i,theta_i);
  p.FaceColor = 'texturemap';
  colormap(pmkmp(56,'Swtth'));
%   colormap(cmocean('thermal',56));
  caxis([-2 12]);
  p.EdgeColor = 'none';         
  alpha(p,0.8);


  % %%% Positive vorticity
  % fv = isosurface(XX,YY,ZZ,vort,0.1);
  % p = patch(fv);
  % p.FaceColor = 'red';
  % p.EdgeColor = 'none';
  % alpha(p,0.5);
  % 
  % %%% Negative vorticity
  % fv = isosurface(XX,YY,ZZ,vort,-0.1);
  % p = patch(fv);
  % p.FaceColor = 'blue';
  % p.EdgeColor = 'none';
  % alpha(p,0.5);

  %%% Isopycnal
  fv = isosurface(XX(:,2:end-1,:),YY(:,2:end-1,:),ZZ(:,2:end-1,:),theta(:,2:end-1,:),theta_plot);
  p = patch(fv);
  % p.FaceColor = 'blue';
  p.FaceColor = [219 110 110]/255;
  p.EdgeColor = 'none';
  alpha(p,0.5);

  save('grfpplot.mat','XX','YY','ZZ','theta','SHELFICEtopo','bathy')

  hold off;
  %pause(1000)

  %%% Decorations
%view(-195,180);
  axis tight;
  xlabel('x (km)','interpreter','latex');
  ylabel('y (km)','interpreter','latex');
  zlabel('z (m)','interpreter','latex');
  %set(gca,'XLim',[-2000 2000]);
  set(gca,'XLim',[-150 150]);
  set(gca,'XTick',[-2000:1000:2000]);
  %set(gca,'YLim',[0 2500]);
  set(gca,'YLim',[0 400]);
  set(gca,'YTick',[0:1000:2500]);
  %set(gca,'ZLim',[-4000 0]);
  set(gca,'ZLim',[-2000 0]);
  set(gca,'ZTick',[-4000:2000:0]);
  set(gca,'FontSize',fontsize);
  pbaspect([Lx/Ly 1 1]);
  set(gca,'Position',[ 0.07   0.07    0.8650    0.8550]);
  handle = colorbar;
  % set(handle,'Position',[0.9337 0.1100 0.02 0.8150]);
  set(handle,'Position',[0.9199    0.6983    0.0201    0.2117]);
  annotation('textbox',[0.01 0.05 0.05 0.01],'String',figlabel,'FontSize',fontsize+2,'LineStyle','None','interpreter','latex');
  annotation('textbox',[0.78 0.9 0.15 0.01],'String',{'Potential';'temperature ($^\circ$C)'},'FontSize',fontsize+2,'LineStyle','None','interpreter','latex');
  title(['t = ',num2str(round(tt),'%.3d'),' days']);
  camlight('headlight');
  lightangle(0,15);
  lighting gouraud;
  
  waitfor(fig);
  M(n) = getframe(gcf);
  
end

%%% Write as .mp4 file
vw = VideoWriter('~/Projects/MITgcm_ISC/pics/theta200','Motion JPEG AVI');
vw.FrameRate=10;
open(vw);
for m = M
  writeVideo(vw,m);
end
close(vw);

