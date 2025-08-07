% Updated script to generate a plot of model bathymetry along with observational data station locations

% created by DJP 10/23

close all
clear

% Load grid file
Grd = ncstruct('./Bering10K_extended_grid.nc');
[nxi, neta] = size(Grd.h);

% Define shelf boundary
isebs = Grd.surveystrata_comboeast < 100;

[x,y] = mask2poly(Grd.lon_psi,Grd.lat_psi,isebs(2:end-1,2:end-1));

BB = shaperead('./BBRKC_Shape_Files/BristolBay.shp');

% Convert Longitude to positive degrees East
BB.BoundingBox(:,1) = BB.BoundingBox(:,1) + 360;

BB.X = BB.X + 360;

% Make masks
[in_BB,on_BB] = inpolygon(Grd.lon_rho,Grd.lat_rho,BB.X,BB.Y);
isBB = in_BB == 1; 
[x_BB,y_BB] = mask2poly(Grd.lon_psi,Grd.lat_psi,isBB(2:end-1,2:end-1));

% Plot figure
figure(1); set(gcf,'Color',[1 1 1]);
m_proj('Mercator','long',[min(min(Grd.lon_rho)) max(max(Grd.lon_rho))],'lat',[min(min(Grd.lat_rho)) max(max(Grd.lat_rho))])
m_pcolor(Grd.lon_rho,Grd.lat_rho,Grd.h);
shading flat
hold on
m_line(x,y,'color','k','Linewidth',3)
m_line(x_BB,y_BB,'color','k','Linewidth',2,'Linestyle','--')
[C200,h200] = m_contour(Grd.lon_rho, Grd.lat_rho, Grd.h, [50 100 200], 'LineWidth', 0.5, 'Color', [0 0 0]);
cb = colorbar;
ylabel(cb,'Depth (m)','Fontsize',14)
caxis([0 8000])
colormap(flipud(cmocean('-deep')))
m_gshhs_i('patch',[0.5 0.5 0.5]);
m_grid('xtick',[160:20:210],'ytick',[50:8:66],'linestyle','none','color','k','Fontsize',6);
title('Bering10K Model Bathymetry','Fontsize',16)
set(gca,'SortMethod','ChildOrder') % To avoid warning message with export_fig


