%%%
%%% plot_model.m
%%%
%%%
function plot_model(expdir,expname,outputname,boundary)

    close all;

    %%% Add path
    %addpath /Users/csi/MITgcm_UC/analysis_uc/functions;
    %addpath /Users/csi/MITgcm_UC/analysis_uc/colormaps;
    %addpath /Users/csi/MITgcm_UC/analysis_uc/colormaps/cmocean/;
    %addpath /Users/csi/Software/eos80_legacy_gamma_n/;
    %addpath /Users/csi/Software/eos80_legacy_gamma_n/library/;
    addpath /home/garrett/Documents/GSW;
    addpath /home/garrett/Documents/GSW/library;
    addpath /home/garrett/Documents/freezecolors;
    %addpath /Users/csi/Software/gsw_matlab_v3_06_11/library/;

    n =1; % Load the reference experiment
%     expdir = EXPDIR{nn};
%     year = YEAR{nn};
    expdir
    loadexp;
    fontsize = 25;
    ncolor=250; % Number of color contours
    m1km = 1000;
    isotherm = true;
    if isotherm
        load_data;
    end


    % bathycolor = hex2rgb('#958872');
    bathycolor = [188 176 153]/255;
    icedraftcolor =  [123 161 210]/255;
    icetopcolor = [123 161 210]/255;
    boxcolor = [0.85 0.85 0.85];
    % isothermcolor = [87 151 246]/255;
    isothermcolor =  [238 131 70]/255;
    
    %%% Calculate CDW layer properties
%     prodir = [expdir expname '/'];
%     load([prodir '/' expname '_tavg_5yrs.mat'], 'THETA','SALT','UVEL','VVEL','VVELTH','ETAN')
%     calc_basics;

    %%% Read snapshot data
%     theta_inst(hFacC==0) = NaN; %%% Remove topography
%     salt_inst(hFacC==0) = NaN; 

    %%% Select potential temperature surface
    theta_plot = 0.5;


    %%% bathymetry and icedraft
    Hicefront = 200; %%% Depth of ice shelf frace
    h=bathy;
    fid = fopen(fullfile(exppath,'input','SHELFICEtopoFile.bin'),'r','b');
    icedraft = fread(fid,[Nx Ny],'real*8');
    
    

    %%% Grid
    [Y,X] = meshgrid(yy/1000,xx/1000);
    [YY,XX,ZZ]=meshgrid(yy/1000,xx/1000,zz/1000);
    [ZZZ,YYY] = meshgrid(zz/1000,yy/1000);

    if isotherm &&  boundary
        idx_1 = 1;
        BC_u = squeeze(uu(idx_1,:,:));
        BC_u(BC_u==0) = NaN;
        BC_t = squeeze(tt(end,:,:));
        BC_s = squeeze(ss(idx_1,:,:));
    end

    %%% Extract zonal boundary values
    idx_1 = 10;
    %%% Calculate the restoring neutral density at the zonal boundaries
    lon_sec = -115;
    lat_sec = -71;
    %%% Plot bathymetry and ice draft
    set(0,'DefaultFigureVisible','off')
    fig = figure(1)
    set(gcf,'Position',[1  107 1475 1139])
    clf;    
    
    %%% Bathymetry  
    hfilled = -h/1000;
    hfilled(:,end)=max(max(abs(ZZZ)));
    p = surface(X,Y,hfilled);
    p.FaceColor = [11*16+9 9*16+12 6*16+11]/255;
    p.FaceColor = bathycolor;
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
    p = surface(X(:,2:end-1),Y(:,2:end-1),-icedraft_plot(:,2:end-1)/1000);
    p.FaceColor = icedraftcolor;
    p.EdgeColor = 'none';
    alpha(p,1);
    p = surface(X(:,2:end-1),Y(:,2:end-1),-icetop_plot(:,2:end-1)/1000);
    p.FaceColor = icetopcolor;
    p.EdgeColor = 'none';
    alpha(p,1);
    if isotherm && boundary
        ttmod = tt;
        ttmod(5,:,:)=NaN;

        ttmod(xx<-50000,:,:)=NaN;
        ttmod(:,:,-zz>1000)=NaN;

        is = isosurface(XX,YY,-ZZ,ttmod,0);
        is= patch(is)

        is.FaceColor = isothermcolor;
        is.EdgeColor = 'none';
        alpha(is,0.5);
        p_bct = surface(xx(end)/1000*ones(size(YYY)),YYY,-ZZZ,BC_t);
        colormap(cmocean('thermal'));
    %     colormap(cmocean('balance'))
        %colormap(cmocean('diff'))
        alpha(p_bct,1);
        %colormap(colormap(cmocean('balance',ncolor)))
        caxis([-2 1]);
        handle_tt = colorbar(gca);
        set(handle_tt,'TickLabels', [-2,0,1 ],'Ticks', [-2,0,1 ]);
        set(handle_tt,'Position',[0.75    0.3    0.0045    0.15]);
        annotation('textbox',[0.6 0.425 0.15 0.01],'String',{'Restoring';'temperature';['(' char(176) 'C)']},'FontSize',fontsize-1,'LineStyle','None','horizontalAlignment','right');
        p_bct.FaceColor = 'texturemap';
        p_bct.EdgeColor = 'none';         
    end
    %freezeColors;

    %%% Decorations
    hold off;
%     view(-219,47);
    view(-223+90,16);
    xlabel('x (km)');
    ylabel('y (km)');
    zlabel('Depth (km)');
    set(gca,'FontSize',fontsize);
    set(gca,'ZTick',[0:0.5:2],'fontsize',fontsize-1);
    set(gca,'YLim',[0 400]);
    set(gca,'YTick',[0:100:400],'fontsize',fontsize-1);
    set(gca,'XLim',[-300 300]);
    set(gca,'XTick',[-300:100:300],'fontsize',fontsize-1);
    set(gca, 'ZDir','reverse')
    axis tight;
    pbaspect([Lx/Ly 1 1]);
    camlight('headlight');
    lightangle(140+90,-34);
    lighting flat;
    box on;
    grid off;

     figdir = '/home/garrett/Projects/HUB/paperfigures/';
    [figdir outputname]
     print('-dpng','-r200',[figdir outputname]);
    
end
    

