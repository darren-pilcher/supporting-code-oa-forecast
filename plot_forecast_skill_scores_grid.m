% Script used to make the grid of ACC values using .mat files calculated in plot_forecast_skill_scores.m
% Created by DJP 01/25

close all
clear

% Load .mat files of ACC stats from forecast
load('./pH_surf_stats.mat');
load('./pH_bot_stats.mat');
load('./arag_surf_stats.mat');
load('./arag_bot_stats.mat');

% Need to append on extra row and column for pcolor grid plotting
pH_surf_stats.ACC(:,10) = pH_surf_stats.ACC(:,9); pH_surf_stats.ACC(3,:) = pH_surf_stats.ACC(2,:);
pH_surf_stats.gt0(:,10) = pH_surf_stats.gt0(:,9); pH_surf_stats.gt0(3,:) = pH_surf_stats.gt0(2,:);
pH_surf_stats.pc5(:,10) = pH_surf_stats.pc5(:,9); pH_surf_stats.pc5(3,:) = pH_surf_stats.pc5(2,:);

pH_bot_stats.ACC(:,10) = pH_bot_stats.ACC(:,9); pH_bot_stats.ACC(3,:) = pH_bot_stats.ACC(2,:);
pH_bot_stats.gt0(:,10) = pH_bot_stats.gt0(:,9); pH_bot_stats.gt0(3,:) = pH_bot_stats.gt0(2,:);
pH_bot_stats.pc5(:,10) = pH_bot_stats.pc5(:,9); pH_bot_stats.pc5(3,:) = pH_bot_stats.pc5(2,:);

arag_surf_stats.ACC(:,10) = arag_surf_stats.ACC(:,9); arag_surf_stats.ACC(3,:) = arag_surf_stats.ACC(2,:);
arag_surf_stats.gt0(:,10) = arag_surf_stats.gt0(:,9); arag_surf_stats.gt0(3,:) = arag_surf_stats.gt0(2,:);
arag_surf_stats.pc5(:,10) = arag_surf_stats.pc5(:,9); arag_surf_stats.pc5(3,:) = arag_surf_stats.pc5(2,:);

arag_bot_stats.ACC(:,10) = arag_bot_stats.ACC(:,9); arag_bot_stats.ACC(3,:) = arag_bot_stats.ACC(2,:);
arag_bot_stats.gt0(:,10) = arag_bot_stats.gt0(:,9); arag_bot_stats.gt0(3,:) = arag_bot_stats.gt0(2,:);
arag_bot_stats.pc5(:,10) = arag_bot_stats.pc5(:,9); arag_bot_stats.pc5(3,:) = arag_bot_stats.pc5(2,:);

% Load .mat files of ACC stats from persistence forecast
load('./Pc_pH_surf_stats.mat');
load('./Pc_pH_bot_stats.mat');
load('./Pc_arag_surf_stats.mat');
load('./Pc_arag_bot_stats.mat');

% Need to append on extra row and column for pcolor grid plotting
Pc_pH_surf_stats.ACC(:,10) = Pc_pH_surf_stats.ACC(:,9); Pc_pH_surf_stats.ACC(3,:) = Pc_pH_surf_stats.ACC(2,:);
Pc_pH_surf_stats.gt0(:,10) = Pc_pH_surf_stats.gt0(:,9); Pc_pH_surf_stats.gt0(3,:) = Pc_pH_surf_stats.gt0(2,:);

Pc_pH_bot_stats.ACC(:,10) = Pc_pH_bot_stats.ACC(:,9); Pc_pH_bot_stats.ACC(3,:) = Pc_pH_bot_stats.ACC(2,:);
Pc_pH_bot_stats.gt0(:,10) = Pc_pH_bot_stats.gt0(:,9); Pc_pH_bot_stats.gt0(3,:) = Pc_pH_bot_stats.gt0(2,:);

Pc_arag_surf_stats.ACC(:,10) = Pc_arag_surf_stats.ACC(:,9); Pc_arag_surf_stats.ACC(3,:) = Pc_arag_surf_stats.ACC(2,:);
Pc_arag_surf_stats.gt0(:,10) = Pc_arag_surf_stats.gt0(:,9); Pc_arag_surf_stats.gt0(3,:) = Pc_arag_surf_stats.gt0(2,:);

Pc_arag_bot_stats.ACC(:,10) = Pc_arag_bot_stats.ACC(:,9); Pc_arag_bot_stats.ACC(3,:) = Pc_arag_bot_stats.ACC(2,:);
Pc_arag_bot_stats.gt0(:,10) = Pc_arag_bot_stats.gt0(:,9); Pc_arag_bot_stats.gt0(3,:) = Pc_arag_bot_stats.gt0(2,:);

cmap = [customcolormap(linspace(0,1,11), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691'})];

% Plot each variable-depth as 1 separate subplot
figure(1); set(gcf, 'units','centimeters','position',[10 40 17 10]); set(gcf,'Color',[1 1 1]);
ax(1) = subplot(4,2,1);
pcolor(pH_surf_stats.ACC)
hold on
tmp = [1.5:1:10.5]; tmp(~pH_surf_stats.gt0(1,:)) = nan;
scatter(tmp,repmat(1.5,[1 10]),60,'k');
tmp = [1.5:1:10.5]; tmp(~pH_surf_stats.gt0(2,:)) = nan;
scatter(tmp,repmat(2.5,[1 10]),60,'k');
tmp = [1.5:1:10.5]; tmp(~pH_surf_stats.pc5(1,:)) = nan;
scatter(tmp,repmat(1.5,[1 10]),60,'k+');
tmp = [1.5:1:10.5]; tmp(~pH_surf_stats.pc5(2,:)) = nan;
scatter(tmp,repmat(2.5,[1 10]),60,'k+');
axis([1 10 1 3])
caxis([0 1])
set(gca,'xtick',[2.5:2:10.5]); set(gca,'xticklabel',[2:2:10])
set(gca,'ytick',[1.5:1:2.5]); set(gca,'yticklabel',{'April' 'May'},'Fontsize',10)
title('Dynamic','Surface pH','Fontsize',10)

ax(2) = subplot(4,2,2);
pcolor(Pc_pH_surf_stats.ACC)
hold on
tmp = [1.5:1:10.5]; tmp(~Pc_pH_surf_stats.gt0(1,:)) = nan;
scatter(tmp,repmat(1.5,[1 10]),60,'k');
tmp = [1.5:1:10.5]; tmp(~Pc_pH_surf_stats.gt0(2,:)) = nan;
scatter(tmp,repmat(2.5,[1 10]),60,'k');
axis([1 10 1 3])
caxis([0 1])
set(gca,'xtick',[2.5:2:10.5]); set(gca,'xticklabel',[2:2:10])
set(gca,'ytick',[1.5:1:2.5]); set(gca,'yticklabel',{'April' 'May'},'Fontsize',10)
title('Persistence','Surface pH','Fontsize',10)

ax(3) = subplot(4,2,3);
pcolor(pH_bot_stats.ACC)
hold on
tmp = [1.5:1:10.5]; tmp(~pH_bot_stats.gt0(1,:)) = nan;
scatter(tmp,repmat(1.5,[1 10]),60,'k');
tmp = [1.5:1:10.5]; tmp(~pH_bot_stats.gt0(2,:)) = nan;
scatter(tmp,repmat(2.5,[1 10]),60,'k');
tmp = [1.5:1:10.5]; tmp(~pH_bot_stats.pc5(1,:)) = nan;
scatter(tmp,repmat(1.5,[1 10]),60,'k+');
tmp = [1.5:1:10.5]; tmp(~pH_bot_stats.pc5(2,:)) = nan;
scatter(tmp,repmat(2.5,[1 10]),60,'k+');
axis([1 10 1 3])
caxis([0 1])
set(gca,'xtick',[2.5:2:10.5]); set(gca,'xticklabel',[2:2:10])
set(gca,'ytick',[1.5:1:2.5]); set(gca,'yticklabel',{'April' 'May'},'Fontsize',10)
title('','Bottom pH','Fontsize',10)

ax(4) = subplot(4,2,4);
pcolor(Pc_pH_bot_stats.ACC)
hold on
tmp = [1.5:1:10.5]; tmp(~Pc_pH_bot_stats.gt0(1,:)) = nan;
scatter(tmp,repmat(1.5,[1 10]),60,'k');
tmp = [1.5:1:10.5]; tmp(~Pc_pH_bot_stats.gt0(2,:)) = nan;
scatter(tmp,repmat(2.5,[1 10]),60,'k');
axis([1 10 1 3])
caxis([0 1])
set(gca,'xtick',[2.5:2:10.5]); set(gca,'xticklabel',[2:2:10])
set(gca,'ytick',[1.5:1:2.5]); set(gca,'yticklabel',{'April' 'May'},'Fontsize',10)
title('','Bottom pH','Fontsize',10)

ax(5) = subplot(4,2,5);
pcolor(arag_surf_stats.ACC)
hold on
tmp = [1.5:1:10.5]; tmp(~arag_surf_stats.gt0(1,:)) = nan;
scatter(tmp,repmat(1.5,[1 10]),60,'k');
tmp = [1.5:1:10.5]; tmp(~arag_surf_stats.gt0(2,:)) = nan;
scatter(tmp,repmat(2.5,[1 10]),60,'k');
tmp = [1.5:1:10.5]; tmp(~arag_surf_stats.pc5(1,:)) = nan;
scatter(tmp,repmat(1.5,[1 10]),60,'k+');
tmp = [1.5:1:10.5]; tmp(~arag_surf_stats.pc5(2,:)) = nan;
scatter(tmp,repmat(2.5,[1 10]),60,'k+');
axis([1 10 1 3])
caxis([0 1])
set(gca,'xtick',[2.5:2:10.5]); set(gca,'xticklabel',[2:2:10])
set(gca,'ytick',[1.5:1:2.5]); set(gca,'yticklabel',{'April' 'May'},'Fontsize',10)
title('','Surface {\Omega}_{arag}','Fontsize',10)

ax(6) = subplot(4,2,6);
pcolor(Pc_arag_surf_stats.ACC)
hold on
tmp = [1.5:1:10.5]; tmp(~Pc_arag_surf_stats.gt0(1,:)) = nan;
scatter(tmp,repmat(1.5,[1 10]),60,'k');
tmp = [1.5:1:10.5]; tmp(~Pc_arag_surf_stats.gt0(2,:)) = nan;
scatter(tmp,repmat(2.5,[1 10]),60,'k');
axis([1 10 1 3])
caxis([0 1])
set(gca,'xtick',[2.5:2:10.5]); set(gca,'xticklabel',[2:2:10])
set(gca,'ytick',[1.5:1:2.5]); set(gca,'yticklabel',{'April' 'May'},'Fontsize',10)
title('','Surface {\Omega}_{arag}','Fontsize',10)

ax(7) = subplot(4,2,7);
pcolor(arag_bot_stats.ACC)
hold on
tmp = [1.5:1:10.5]; tmp(~arag_bot_stats.gt0(1,:)) = nan;
scatter(tmp,repmat(1.5,[1 10]),60,'k');
tmp = [1.5:1:10.5]; tmp(~arag_bot_stats.gt0(2,:)) = nan;
scatter(tmp,repmat(2.5,[1 10]),60,'k');
tmp = [1.5:1:10.5]; tmp(~arag_bot_stats.pc5(1,:)) = nan;
scatter(tmp,repmat(1.5,[1 10]),60,'k+');
tmp = [1.5:1:10.5]; tmp(~arag_bot_stats.pc5(2,:)) = nan;
scatter(tmp,repmat(2.5,[1 10]),60,'k+');
axis([1 10 1 3])
caxis([0 1])
set(gca,'xtick',[2.5:2:10.5]); set(gca,'xticklabel',[2:2:10])
set(gca,'ytick',[1.5:1:2.5]); set(gca,'yticklabel',{'April' 'May'},'Fontsize',10)
xlabel('Lead time (months)','Fontsize',10)
title('','Bottom {\Omega}_{arag}','Fontsize',10)

ax(8) = subplot(4,2,8);
pcolor(Pc_arag_bot_stats.ACC)
hold on
tmp = [1.5:1:10.5]; tmp(~Pc_arag_bot_stats.gt0(1,:)) = nan;
scatter(tmp,repmat(1.5,[1 10]),60,'k');
tmp = [1.5:1:10.5]; tmp(~Pc_arag_bot_stats.gt0(2,:)) = nan;
scatter(tmp,repmat(2.5,[1 10]),60,'k');
axis([1 10 1 3])
caxis([0 1])
set(gca,'xtick',[2.5:2:10.5]); set(gca,'xticklabel',[2:2:10])
set(gca,'ytick',[1.5:1:2.5]); set(gca,'yticklabel',{'April' 'May'},'Fontsize',10)
xlabel('Lead time (months)','Fontsiz',10)
title('','Bottom {\Omega}_{arag}','Fontsize',10)
colormap(cmap)

cb = colorbar; set(cb,'position',[.915 ax(8).Position(2) .025 .77],'orientation','vertical','Fontsize',8)
ylabel(cb,'ACC','Fontsize',10)








