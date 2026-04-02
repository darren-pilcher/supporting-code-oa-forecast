% Script used to plot forecast skill scores for temp,salt,DIC, and TA for underyling mechanisms of predictability

% Created by DJP 03/26

close all
clear

DIR = './';

% Load skill metrics calcualted and saved in plot_forecast_skill_scores.m
VAR = {'temp', 'salt','alkalinity','TIC','arag','pH'};
for v = 1:length(VAR)
	load([DIR,char(VAR(v)),'_surf_stats.mat']);
    load([DIR,char(VAR(v)),'_bot_stats.mat']);
	load([DIR,'Pc_',char(VAR(v)),'_surf_stats.mat']);
    load([DIR,'Pc_',char(VAR(v)),'_bot_stats.mat']);
end

% Calculate correlations between mech vars and pH, arag 
for s = 1:2
	arag_corr.surf(s,1) = corr(arag_surf_stats.ACC(s,:)',temp_surf_stats.ACC(s,:)');
    arag_corr.surf(s,2) = corr(arag_surf_stats.ACC(s,:)',salt_surf_stats.ACC(s,:)');
	arag_corr.surf(s,3) = corr(arag_surf_stats.ACC(s,:)',TIC_surf_stats.ACC(s,:)');
	arag_corr.surf(s,4) = corr(arag_surf_stats.ACC(s,:)',alkalinity_surf_stats.ACC(s,:)');

    arag_corr.bot(s,1) = corr(arag_bot_stats.ACC(s,:)',temp_bot_stats.ACC(s,:)');
    arag_corr.bot(s,2) = corr(arag_bot_stats.ACC(s,:)',salt_bot_stats.ACC(s,:)');
    arag_corr.bot(s,3) = corr(arag_bot_stats.ACC(s,:)',TIC_bot_stats.ACC(s,:)');
    arag_corr.bot(s,4) = corr(arag_bot_stats.ACC(s,:)',alkalinity_bot_stats.ACC(s,:)');

    pH_corr.surf(s,1) = corr(pH_surf_stats.ACC(s,:)',temp_surf_stats.ACC(s,:)');
    pH_corr.surf(s,2) = corr(pH_surf_stats.ACC(s,:)',salt_surf_stats.ACC(s,:)');
    pH_corr.surf(s,3) = corr(pH_surf_stats.ACC(s,:)',TIC_surf_stats.ACC(s,:)');
    pH_corr.surf(s,4) = corr(pH_surf_stats.ACC(s,:)',alkalinity_surf_stats.ACC(s,:)');

    pH_corr.bot(s,1) = corr(pH_bot_stats.ACC(s,:)',temp_bot_stats.ACC(s,:)');
    pH_corr.bot(s,2) = corr(pH_bot_stats.ACC(s,:)',salt_bot_stats.ACC(s,:)');
    pH_corr.bot(s,3) = corr(pH_bot_stats.ACC(s,:)',TIC_bot_stats.ACC(s,:)');
    pH_corr.bot(s,4) = corr(pH_bot_stats.ACC(s,:)',alkalinity_bot_stats.ACC(s,:)');

end

% Plot figure
figure(1); set(gcf, 'units','centimeters','position',[10 40 29 15]); set(gcf,'Color',[1 1 1]);
titles = ["April Initialized Surface" "May Initialized Surface" "April Initialized Bottom" "May Initialized Bottom"];
for s = 1:2
	subplot(2,2,s)
	ax1 = plot(temp_surf_stats.ACC(s,:),'Linewidth',2);
	hold on
	tmp = temp_surf_stats.ACC(s,:); tmp(~temp_surf_stats.pc5(s,:)) = nan;
	plot(tmp,'k','Marker','d','MarkerFaceColor',ax1.Color,'LineStyle','none')
	ax2 = plot(salt_surf_stats.ACC(s,:),'Linewidth',2);
	tmp = salt_surf_stats.ACC(s,:); tmp(~salt_surf_stats.pc5(s,:)) = nan;
	plot(tmp,'k','Marker','d','MarkerFaceColor',ax2.Color,'LineStyle','none')
	ax3 = plot(TIC_surf_stats.ACC(s,:),'Linewidth',2);
	tmp = TIC_surf_stats.ACC(s,:); tmp(~TIC_surf_stats.pc5(s,:)) = nan;
	plot(tmp,'k','Marker','d','MarkerFaceColor',ax3.Color,'LineStyle','none')
	ax4 = plot(alkalinity_surf_stats.ACC(s,:),'Color',[0.4940 0.1840 0.5560],'Linewidth',2);
	tmp = alkalinity_surf_stats.ACC(s,:); tmp(~alkalinity_surf_stats.pc5(s,:)) = nan;
	plot(tmp,'k','Marker','d','MarkerFaceColor',ax4.Color,'LineStyle','none')
	plot(ones(1,9)*0.5,'k--','Linewidth',1)
	axis([1 9 -.05 1])
	ylabel('ACC','Fontsize',14)
	title(titles(s),'Fontsize',16)

	subplot(2,2,s+2)
	ax1 = plot(temp_bot_stats.ACC(s,:),'Linewidth',2);
	hold on
	tmp = temp_bot_stats.ACC(s,:); tmp(~temp_bot_stats.pc5(s,:)) = nan;
	plot(tmp,'k','Marker','d','MarkerFaceColor',ax1.Color,'LineStyle','none')
	ax2 = plot(salt_bot_stats.ACC(s,:),'Linewidth',2);
	tmp = salt_bot_stats.ACC(s,:); tmp(~salt_bot_stats.pc5(s,:)) = nan;
	plot(tmp,'k','Marker','d','MarkerFaceColor',ax2.Color,'LineStyle','none')
	ax3 = plot(TIC_bot_stats.ACC(s,:),'Linewidth',2);
	tmp = TIC_bot_stats.ACC(s,:); tmp(~TIC_bot_stats.pc5(s,:)) = nan;
	plot(tmp,'k','Marker','d','MarkerFaceColor',ax3.Color,'LineStyle','none')
	ax4 = plot(alkalinity_bot_stats.ACC(s,:),'Color',[0.4940 0.1840 0.5560],'Linewidth',2);
	tmp = alkalinity_bot_stats.ACC(s,:); tmp(~alkalinity_bot_stats.pc5(s,:)) = nan;
	plot(tmp,'k','Marker','d','MarkerFaceColor',ax4.Color,'LineStyle','none')
	plot(ones(1,9)*0.5,'k--','Linewidth',1)
	axis([1 9 -.05 1])
	ylabel('ACC','Fontsize',14)
        xlabel('Lead time (months)','Fontsize',14)
	title(titles(s+2),'Fontsize',16)
end
legend([ax1,ax2,ax3,ax4],'Temp','Salt','DIC','TA','Location','Southwest')


