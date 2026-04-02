% Script used to take the spatial ACC values for 2 variables and run a spatial correlation

% Created by DJP 03/26

close all
clear

DIR = './';

% Select which variables to compare
VAR1 = 'arag';
VAR2 = 'temp';

% Load grid file
Grd = ncstruct('./Bering10K_extended_grid.nc');

% Load .mat files calculated using plot_forecast_skill_scores_spatial.m
load([DIR,VAR1,'_surf_spatial_ACC_may.mat']); var1_surf = eval([VAR1,'_surf_spatial_ACC_may']);
load([DIR,VAR1,'_bot_spatial_ACC_may.mat']); var1_bot = eval([VAR1,'_bot_spatial_ACC_may']);

load([DIR,VAR2,'_surf_spatial_ACC_may.mat']); var2_surf = eval([VAR2,'_surf_spatial_ACC_may']);
load([DIR,VAR2,'_bot_spatial_ACC_may.mat']); var2_bot = eval([VAR2,'_bot_spatial_ACC_may']);
clear *_ACC*

% Calculate spatial correltion
for y = 1:size(var1_surf,2)
	for x = 1:size(var1_surf,1)
		vars_surf_corr(x,y) = corr(squeeze(var1_surf(x,y,:)),squeeze(var2_surf(x,y,:)));
                vars_bot_corr(x,y) = corr(squeeze(var1_bot(x,y,:)),squeeze(var2_bot(x,y,:)));

	end
end

% Spatial masks
isebs = Grd.surveystrata_updated < 100 & isnan(vars_surf_corr)==0;

[X,Y] = mask2poly(Grd.lon_psi,Grd.lat_psi,isebs(2:end-1,2:end-1));

cmap = customcolormap(linspace(0,1,11), {'#523107','#523107','#bf812c','#e2c17e','#f3e9c4','#f6f4f4','#cae9e3','#81cdc1','#379692','#01665e','#003d2e'}); % brown-white-pool

vars_surf_corr_ebs = local(vars_surf_corr,isebs,'weight',Grd.area_feast,'omitnan'); 
vars_bot_corr_ebs = local(vars_bot_corr,isebs,'weight',Grd.area_feast,'omitnan');

% Plot figures
figure(1); set(gcf, 'units','centimeters','position',[10 20 25 13]); set(gcf,'Color',[1 1 1]);
ax1 = subplot(1,2,1);
m_proj('Mercator','long',[177 205],'lat',[53 66])
fi = find(Grd.h > 1500); tmp_plot = vars_surf_corr; tmp_plot(fi) = nan;
[C,h] = m_contourf(Grd.lon_rho,Grd.lat_rho,tmp_plot,[-1:.025:1]);
set(h,'LineColor','none')
hold on
m_line(X,Y,'color','k','Linewidth',2)
[C200,h200] = m_contour(Grd.lon_rho, Grd.lat_rho, Grd.h, [50 100 200], 'LineWidth', 0.5, 'Color', [0 0 0]);
caxis([-1 1])
m_gshhs_i('patch',[0.5 0.5 0.5]);
m_grid('xtick',[160:20:210],'ytick',[50:8:66],'linestyle','none','color','k','Fontsize',6);
title(['Surface ',VAR1,' ',VAR2,' Correlation (Avg. = ',num2str(vars_surf_corr_ebs,'%#.2f'),')'],'Fontsize',12)

ax2 = subplot(1,2,2);
m_proj('Mercator','long',[177 205],'lat',[53 66])
fi = find(Grd.h > 1500); tmp_plot = vars_bot_corr; tmp_plot(fi) = nan;
[C,h] = m_contourf(Grd.lon_rho,Grd.lat_rho,tmp_plot,[-1:.025:1]);
set(h,'LineColor','none')
hold on
m_line(X,Y,'color','k','Linewidth',2)
[C200,h200] = m_contour(Grd.lon_rho, Grd.lat_rho, Grd.h, [50 100 200], 'LineWidth', 0.5, 'Color', [0 0 0]);
caxis([-1 1])
m_gshhs_i('patch',[0.5 0.5 0.5]);
m_grid('xtick',[160:20:210],'ytick',[50:8:66],'linestyle','none','color','k','Fontsize',6);
title(['Bottom ',VAR1,' ',VAR2,' Correlation (Avg. = ',num2str(vars_bot_corr_ebs,'%#.2f'),')'],'Fontsize',12)

colormap(flipud(cmap)); 
cb = colorbar;
set(cb,'position',[ax1.Position(1) ax2.Position(2)+.01 ax2.Position(1)+ax2.Position(3)-ax1.Position(1) .05],'orientation','horizontal')
xlabel(cb,'Correlation Coefficient','Fontsize',12)



