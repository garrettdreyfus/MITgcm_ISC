% %%%
%%% anim.m
%%%
%%% Reads diagnostic output from MITgcm and makes a movie of it. Requires
%%% that the experiment be generated using the 'newexp' package.
%%%

%%% NOTE: Doesn't account for u/v gridpoint locations

%%% Set true if plotting on a Mac
mac_plots = 0;

%%% Read experiment data
loadexp;

%%% Select diagnostic variable to animate
diagnum = 3;
outfname = diag_fileNames{1,diagnum};

%%% Data index in the output data files
outfidx = 1;

%%% If set true, plots a top-down view of the field in a given layer.
%%% Otherwise plots a side-on view of the zonally-averaged field
xyplot = 0;

%%% Vertical layer index to use for top-down plots
xylayer = 1;

%%% Set true to plot the field in the lowest active cell at each horizontal
%%% location
botplot = 0;

%%% Set true to plot the field in the topmost wet cell at each horizontal
%%% location
topplot = 0;

%%% Set true to plot the field in the middle of the water column at each
%%% horizontal location
midplot = 0;

%%% Set true for a zonal average
yzavg = 1;

%%% Layer to plot in the y/z plane
yzlayer = 80;

%%% Specify color range
set_crange = 1;


% crange = [-2.2 -1.6]; %/%% Filchner temp
crange = [-2 1]; %%%temp
% crange = [34.1 34.7]; %%% salinity
% crange = [33.9 34.3]; %%% surface salinity
% crange = [0 10]; %%%% for KPP hbl
% crange = [0 1]; %%% For sea ice area
% crange = [-.3 .3]; %%% For velocities or stresses
% crange = [-1 1]*1e-4; %%% For freshwater fluxes
% crange =[-100 100]; %%% Qnet
% crange = [-300 300]; %%% swnet
% crange = [0 1]; %%% SI thickness
% crange = [-0.01 0.01];

% cmap = pmkmp(100,'Swtth');
% cmap = cmocean('haline',100);
cmap = cmocean('thermal',100);
% cmap = cmocean('ice',100);
% cmap = haxby;
% cmap = jet(200);
% cmap = redblue(100);

% titlestr = 'Bottom salinity (g/kg)';
% titlestr = 'Sea ice concentration';
titlestr = '';

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(diagnum));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((0:nDumps)*dumpFreq/deltaT);
% dumpIters = dumpIters(dumpIters > nIter0);

%%% Map of shallowest and deepest wet grid cells
kmax = ones(Nx,Ny);
kmin = ones(Nx,Ny);
for i=1:Nx
  for j=1:Ny
    idx = find(squeeze(hFacC(i,j,:))>0);
    if (~isempty(idx))
      kmin(i,j) = min(idx);
      kmax(i,j) = max(idx);
    end
  end
end

%%% Mesh grids for plotting
if (xyplot)

  [YY,XX] = meshgrid(yy,xx);
  if (botplot || topplot || midplot)      
    kn = ones(Nx,Ny);
    kp= ones(Nx,Ny);
    wn = 0.5*ones(Nx,Ny);
    wp = 0.5*ones(Nx,Ny);
    for i=1:Nx
      for j=1:Ny
        if (sum(hFacC(i,j,:),3)==0)
          continue;
        end
        zmid = 0.5 * (SHELFICEtopo(i,j) + bathy(i,j));
        kmid = max(find(squeeze(zz)>zmid));
        if (isempty(kmid))
          continue;
        end
        kp(i,j) = kmid;
        kn(i,j) = kp(i,j) + 1;
        wp(i,j) = (zmid-zz(kn(i,j))) / (zz(kp(i,j))-zz(kn(i,j)));
        wn(i,j) = 1 - wp(i,j);
      end
    end
  end

else
  
  %%% Create mesh grid with vertical positions adjusted to sit on the bottom
  %%% topography and at the surface
  [ZZ,YY] = meshgrid(zz,yy);  
%   for j=1:Ny
%     if (yzavg)
%       hFacC_col = squeeze(hFacC(:,j,:));    
%       hFacC_col = max(hFacC_col,[],1);    
%     else
%       hFacC_col = squeeze(hFacC(yzlayer,j,:))';
%     end
%     zz_topface = zz(kmin(i,j))-(0.5-hFacC_col(kmin(i,j)))*delR(kmin(i,j));
%     zz_botface = zz(kmax(i,j))+(0.5-hFacC_col(kmax(i,j)))*delR(kmax(i,j));
%     ZZ(j,kmin(i,j)) = zz_topface;
%     if (kmax(i,j)>1)
%       ZZ(j,kmax(i,j)) = zz_botface;
%     end
%   end
%   
end

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 26;
if (mac_plots)  
  framepos = [0 scrsz(4)/2 scrsz(3)/1.3 scrsz(4)];
  plotloc = [0.17 0.3 0.62 0.7];
else
  plotloc = [0.11 0.15 0.76 0.73];
  framepos = [327    80   941   885];
end

%%% Set up the figure
handle = figure(10);
% set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
set(gcf,'color','w');
M = moviein(nDumps);

Amean = [];
tyears = [];
Aprev = [];
A = [];
% for n=1
for n=1:length(dumpIters)
% for n=40:length(dumpIters)
% for n=length(dumpIters)-1

  t = dumpIters(n)*deltaT/86400/365;
  
  if (~isempty(A))
    Aprev = A;
  end
  
  A = rdmdsWrapper(fullfile(exppath,'results',outfname),dumpIters(n));          
  if (isempty(A))
    error(['Ran out of data at t=,',num2str(tyears(n)),' years']);
  end  
  
  if (isempty(Aprev))
    Aprev = 0*A;
  end
  
%   if (n > 1)
%     ['diff = ',num2str(max(max(max(abs(A-Aprev)))))]
%   end
 
  tyears(n) = t;
  

  
  
%   DX = repmat(delX',[1 Ny Nr]);
%   DY = repmat(delY,[Nx 1 Nr]);
%   DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);  
%   Amean(n) = sum(sum(sum(A.*DX.*DY.*DZ.*hFacC)))/sum(sum(sum(DX.*DY.*DZ.*hFacC)));
  
%   A(hFacC==0)   = NaN;
[min(A(:),[],'omitnan') max(A(:),[],'omitnan')]
  
    
%   Axy = sum(A.*DY.*DZ.*hFacW,3);
  
  %%% x/y plot
  if (xyplot)
    
    if (botplot)      
      FF = zeros(Nx,Ny);
      for i=1:Nx
        for j=1:Ny
          FF(i,j) = A(i,j,kmax(i,j));
        end
      end
      
    elseif (topplot)
      FF = zeros(Nx,Ny);
      for i=1:Nx
        for j=1:Ny
          FF(i,j) = A(i,j,kmin(i,j));
        end

      end
    elseif (midplot)
      FF = zeros(Nx,Ny);
      for i=1:Nx
        for j=1:Ny
          FF(i,j) = wp(i,j)*A(i,j,kp(i,j)) + wn(i,j)*A(i,j,kn(i,j));
        end
      end

    else
      FF = squeeze(A(:,:,xylayer,outfidx));        
      FF(hFacC(:,:,xylayer)==0) = NaN;
    end
    
    FF(sum(hFacC,3)==0) = NaN;
    
%     contourf(XX,YY,FF,100,'EdgeColor','None');  
    pcolor(XX,YY,FF);
    shading interp;
    xlabel('x (km)');
    ylabel('y (km)');
    
    
    
    
  %%% y/z zonally-averaged plot
  else
    
    if (yzavg)
      Ayz = squeeze(A(:,:,:,outfidx));
      Ayz(hFacC==0) = NaN;
Ayz = squeeze(mean(Ayz,1,'omitnan'));    
    else
      Ayz = squeeze(A(yzlayer,:,:,outfidx));
      Ayz(squeeze(hFacC(yzlayer,:,:))==0) = NaN;
    end
    
    jrange = 1:Ny;
%     [C h] = contourf(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),200,'EdgeColor','None');              
    pcolor(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:));
    shading interp;
    hold on;
%     [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),10,'EdgeColor','k');
%     [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),[-2:0.5:12],'EdgeColor','k');
%     [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),[-0.1:0.005:0.1],'EdgeColor','k');
    if (yzavg)
      h = plot(yy/1000,min(bathy,[],1)/1000,'k','LineWidth',3);  
    else
      h = plot(yy/1000,bathy(yzlayer,:)/1000,'k','LineWidth',3);  
    end
    hold off;
    xlabel('Offshore $y$ (km)','interpreter','latex');
    ylabel('Height $z$ (km)','interpreter','latex');  
  
  end
    
  
  
  %%% Finish the plot
  handle=colorbar;
  caxis(crange);
  colormap(cmap);
  set(handle,'FontSize',fontsize);
  set(gca,'FontSize',fontsize);
  title([titlestr,' at $t=',num2str(tyears(n),'%.1f'),'$ years'],'interpreter','latex');  
  set(gca,'Position',plotloc);
%   annotation('textbox',[0.85 0.05 0.25 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
  M(n) = getframe(gcf);
  
end
