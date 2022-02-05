
  
  
%%% TODO
  use_trough = false;
  H = 4000;
  m1km = 1000;
  Lx = 400*m1km;
  Ly = 400*m1km;
  Nx = 400;
  Ny = 400;
  Nr = 70;
  
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
  
   %%% Zonal grid
  dx = Lx/Nx;  
  xx = (1:Nx)*dx;  
  
  %%% Uniform meridional grid   
  dy = (Ly/Ny)*ones(1,Ny);  
  yy = cumsum((dy + [0 dy(1:end-1)])/2);
 
  %%% Plotting mesh
  [Y,X] = meshgrid(yy,xx);
  
  %%% TODO refine for ice shelf cavity
  %%% Grid spacing increases with depth, but spacings exactly sum to H
  zidx = 1:Nr;
  gamma = 10;  
  alpha = 10;
  dz1 = 2*H/Nr/(alpha+1);
  dz2 = alpha*dz1;
  dz = dz1 + ((dz2-dz1)/2)*(1+tanh((zidx-((Nr+1)/2))/gamma));
  zz = -cumsum((dz+[0 dz(1:end-1)])/2);
  
  
  
  
  
  
  
  
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
  
  

  
  
  
  
  %%% Construct ice she
  icedraft = zeros(Nx,Ny);
  iceidx = find(Y<=Yicefront);  
  icedraft(iceidx) = -Hicefront - (Y(iceidx)-Yicefront)/Yicefront * Hice;
  icedraft(icedraft<h) = h(icedraft<h);
  
  
  
  
  
  
  
  
  fontsize = 12;
  fignum = 1;
      figure(fignum);
    fignum = fignum + 1;
    clf;
    surf(X/1000,Y/1000,h,'EdgeColor','None');   
    xlabel('x (km)');
    ylabel('y (km)');
    zlabel('hb','Rotation',0);
%     plot(Y(1,:),h(1,:));
    title('Model bathymetry');
    set(gca,'fontsize',fontsize+2);
    PLOT = gcf;
    PLOT.Position = [248 284 655 442];  
    
figure(fignum);
fignum = fignum + 1;
clf;    

p = surface(X(:,2:end-1)/1000,Y(:,2:end-1)/1000,h(:,2:end-1));
p.FaceColor = [11*16+9 9*16+12 6*16+11]/255;
p.EdgeColor = 'none';

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