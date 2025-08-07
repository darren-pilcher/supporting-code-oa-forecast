% Script used to calculate the ACC scores for the seasonal forecasts

% General outline here it to 
% 1. Load hindcast output into monthly averaged values for each year to make climatology and persistence forecast
% 2. Load forecast output and also generate monthly averaged values
% 3. Compute anomalies between forecast and hindcast, and then persistence and hindcast


close all
clear

% Select variable
VAR = 'pH';
%VAR = 'arag';

% Select depth
%depth = 'surface5m';
depth = 'bottom5m';

% Detrend or not
detrend = false;

years = [1982:2010];

ens = ["0327","0426"];
mem = ["06", "12", "18"];

% Directories
RUN_DIR_HC = './B10K-K20P19';
RUN_DIR_FC = './CFSR_BCOR/';

% Load grid file
Grd = ncstruct('./Bering10K_extended_grid.nc');

% Load hindcast
F = dir(fullfile(RUN_DIR_HC, ['*average_',VAR,'_',depth,'.nc']));
fname = fullfile({F.folder}, {F.name});

tmp = ncstruct(fname, VAR);
var_tmp = tmp.(VAR);
var_t = ncdateread(fname, 'ocean_time'); clear tmp

% Grab points over selected timeframe
ind_tmp = isbetween(var_t,datetime(years(1),01,01),datetime(years(end),12,31,23,59,59));

Hc.var = var_tmp(:,:,ind_tmp);
Hc.t = var_t(ind_tmp);
clear var*

% Calculate monthly averages in prep for montly climatology
end_day = [31 28 31 30 31 30 31 31 30 31 30 31];

% Loop through each month
for m = 1:12
        day_st = datetime(years,m,1); day_en = datetime(years,m,end_day(m),23,59,59);
        for t = 1:length(years)
                tmp_ind = isbetween(Hc.t,day_st(t),day_en(t));
                Hc.var_mon(:,:,m,t) = squeeze(mean(Hc.var(:,:,tmp_ind),3));
        end
end

% Load forecasts
% Load bottom pH and arag output from forecasts (leave out 1982-1983 for now since they have missing ensemble members)
for t = 1:length(years)-2
        for i = 1:2
                for e = 1:3
                        var_tmp = ncread(sprintf('%sCFSR_BCOR_%d%s%s_average_%s_%s.nc',RUN_DIR_FC,years(t+2),ens(i),mem(e),VAR,depth),VAR);
                        tmp_t = ncdateread(sprintf('%sCFSR_BCOR_%d%s%s_average_%s_%s.nc',RUN_DIR_FC,years(t+2),ens(i),mem(e),VAR,depth),'ocean_time');

			% Monthly averaging
			if i==1
				for m = 1:9
					day_st = datetime(years(t+2),m+3,1); day_en = datetime(years(t+2),m+3,end_day(m+3),23,59,59); 
					tmp_ind = isbetween(tmp_t,day_st,day_en);
					Fc.var_mon_apr(:,:,e,m,t+2) = squeeze(mean(var_tmp(:,:,tmp_ind),3));
				end	
			elseif i==2
				for m = 1:8
					day_st = datetime(years(t+2),m+4,1); day_en = datetime(years(t+2),m+4,end_day(m+4),23,59,59);
                                        tmp_ind = isbetween(tmp_t,day_st,day_en);
                                        Fc.var_mon_may(:,:,e,m,t+2) = squeeze(mean(var_tmp(:,:,tmp_ind),3));
				end	
				% a bit clunky but accounts for Jan of following year
                                years_end = [1982:2011];
                                day_st = datetime(years_end(t+3),1,1); day_en = datetime(years_end(t+3),1,end_day(1),23,59,59);
                                tmp_ind = isbetween(tmp_t,day_st,day_en);
                                Fc.var_mon_may(:,:,e,9,t+2) = squeeze(mean(var_tmp(:,:,tmp_ind),3));
			end
		end
	end
end

% 1982 and 1983 missing (1983032712 1982032718)
% Load May forecasts like usual 
for t = 1:2
	for e = 1:3
		var_tmp = ncread(sprintf('%sCFSR_BCOR_%d%s%s_average_%s_%s.nc',RUN_DIR_FC,years(t),ens(2),mem(e),VAR,depth),VAR);
                tmp_t = ncdateread(sprintf('%sCFSR_BCOR_%d%s%s_average_%s_%s.nc',RUN_DIR_FC,years(t),ens(2),mem(e),VAR,depth),'ocean_time');

		% Monthly averaging
		for m = 1:8
			day_st = datetime(years(t),m+4,1); day_en = datetime(years(t),m+4,end_day(m+4),23,59,59);
                        tmp_ind = isbetween(tmp_t,day_st,day_en);
                        Fc.var_mon_may(:,:,e,m,t) = squeeze(mean(var_tmp(:,:,tmp_ind),3));	
		end
		day_st = datetime(years(t+1),1,1); day_en = datetime(years(t+1),1,end_day(1),23,59,59);
                tmp_ind = isbetween(tmp_t,day_st,day_en);
                Fc.var_mon_may(:,:,e,9,t) = squeeze(mean(var_tmp(:,:,tmp_ind),3));
	end
end

% Load April forecasts
for e = 1:2
	t = 1;
	var_tmp = ncread(sprintf('%sCFSR_BCOR_%d%s%s_average_%s_%s.nc',RUN_DIR_FC,years(t),ens(1),mem(e),VAR,depth),VAR);
        tmp_t = ncdateread(sprintf('%sCFSR_BCOR_%d%s%s_average_%s_%s.nc',RUN_DIR_FC,years(t),ens(1),mem(e),VAR,depth),'ocean_time');

	% Monthly averaging
	for m = 1:9
		day_st = datetime(years(t),m+3,1); day_en = datetime(years(t),m+3,end_day(m+3),23,59,59);
                tmp_ind = isbetween(tmp_t,day_st,day_en);
                Fc.var_mon_apr(:,:,e,m,t) = squeeze(mean(var_tmp(:,:,tmp_ind),3));
        end
end

for e = [1 3]
	t = 2;
        var_tmp = ncread(sprintf('%sCFSR_BCOR_%d%s%s_average_%s_%s.nc',RUN_DIR_FC,years(t),ens(1),mem(e),VAR,depth),VAR);
        tmp_t = ncdateread(sprintf('%sCFSR_BCOR_%d%s%s_average_%s_%s.nc',RUN_DIR_FC,years(t),ens(1),mem(e),VAR,depth),'ocean_time');

        % Monthly averaging
        for m = 1:9
                day_st = datetime(years(t),m+3,1); day_en = datetime(years(t),m+3,end_day(m+3),23,59,59);
                tmp_ind = isbetween(tmp_t,day_st,day_en);
                Fc.var_mon_apr(:,:,e,m,t) = squeeze(mean(var_tmp(:,:,tmp_ind),3));
        end
end


% Fill in nans for missing ensemble members
Fc.var_mon_apr(:,:,3,:,1) = nan; Fc.var_mon_apr(:,:,2,:,2) = nan;

% Detrend hindcast.  Start by calculating trend over annual average values
if (detrend)
	tmp = squeeze(nanmean(Hc.var_mon,3));
	Hc.lin_trend = trend(tmp); clear tmp

% Now remove annual lin trend from monthly hindcast output
	for t = 1:length(years)
		for m = 1:12
			Hc.var_mon(:,:,m,t) = Hc.var_mon(:,:,m,t) - (Hc.lin_trend*t);
		end
	end

% Now detrend forecasts
	for t = 1:length(years)
		for e = 1:3
			for m = 1:9
				Fc.var_mon_apr(:,:,e,m,t) = Fc.var_mon_apr(:,:,e,m,t) - (Hc.lin_trend*t);
			end

			for m = 1:9
				Fc.var_mon_may(:,:,e,m,t) = Fc.var_mon_may(:,:,e,m,t) - (Hc.lin_trend*t);
			end
		end
	end
			

end % detrend

% Compute hindcast monthly anomalies
for t = 1:length(years)
	for m = 1:12
		Hc.var_mon_anom(:,:,m,t) = squeeze(Hc.var_mon(:,:,m,t)) - squeeze(nanmean(Hc.var_mon(:,:,m,:),4));
	end
end

% Compute forecast monthly anomalies
% April
for t = 1:length(years)
	for e = 1:3
		for m = 1:9
			Fc.var_mon_apr_anom(:,:,e,m,t) = squeeze(Fc.var_mon_apr(:,:,e,m,t)) - squeeze(nanmean(Fc.var_mon_apr(:,:,e,m,:),5));
		end
	end
% Also calculate ensemble mean anomalies
	for m = 1:9
		Fc.var_mon_apr_ensmean_anom(:,:,m,t) = squeeze(nanmean(Fc.var_mon_apr(:,:,:,m,t),3)) - squeeze(nanmean(nanmean(Fc.var_mon_apr(:,:,:,m,:),3),5));; 
	end

end

% May
for t = 1:length(years)
        for e = 1:3
                for m = 1:9
                        Fc.var_mon_may_anom(:,:,e,m,t) = squeeze(Fc.var_mon_may(:,:,e,m,t)) - squeeze(nanmean(Fc.var_mon_may(:,:,e,m,:),5));
                end
        end
% Also calculate ensemble mean anomalies
        for m = 1:9
                Fc.var_mon_may_ensmean_anom(:,:,m,t) = squeeze(nanmean(Fc.var_mon_may(:,:,:,m,t),3)) - squeeze(nanmean(nanmean(Fc.var_mon_may(:,:,:,m,:),3),5));;
        end
end

% Compute ACC between forecasts and hidncast
for e = 1:3
	for y = 1:size(Fc.var_mon_apr_anom,2)
		for x = 1:size(Fc.var_mon_apr_anom,1)
			tmp_hc = squeeze(Hc.var_mon_anom(x,y,4:12,:))'; tmp_fc = squeeze(Fc.var_mon_apr_anom(x,y,e,:,:))';
			tmp_corr = corr(tmp_hc,tmp_fc,'rows','complete');
			ACC.Fc_apr(x,y,e,:) = tmp_corr(sub2ind(size(tmp_corr),[1:9],[1:9]));		

                        tmp_hc = squeeze(Hc.var_mon_anom(x,y,5:12,:))'; tmp_fc = squeeze(Fc.var_mon_may_anom(x,y,e,:,:))';
                        tmp_hc_end = [squeeze(Hc.var_mon_anom(x,y,1,2:end));nan];
			tmp_hc = [tmp_hc tmp_hc_end];
			tmp_corr = corr(tmp_hc,tmp_fc,'rows','complete');
                        ACC.Fc_may(x,y,e,:) = tmp_corr(sub2ind(size(tmp_corr),[1:9],[1:9]));

		end
	end
end		

% Same as above but for ensemble mean forecast
for y = 1:size(Fc.var_mon_apr_anom,2)
	for x = 1:size(Fc.var_mon_apr_anom,1)
        	tmp_hc = squeeze(Hc.var_mon_anom(x,y,4:12,:))'; tmp_fc = squeeze(Fc.var_mon_apr_ensmean_anom(x,y,:,:))';
                tmp_corr = corr(tmp_hc,tmp_fc,'rows','complete');
                ACC.Fc_apr_ensmean(x,y,:) = tmp_corr(sub2ind(size(tmp_corr),[1:9],[1:9]));

                tmp_hc = squeeze(Hc.var_mon_anom(x,y,5:12,:))'; tmp_fc = squeeze(Fc.var_mon_may_ensmean_anom(x,y,:,:))';
                tmp_hc_end = [squeeze(Hc.var_mon_anom(x,y,1,2:end));nan];
                tmp_hc = [tmp_hc tmp_hc_end];
		tmp_corr = corr(tmp_hc,tmp_fc,'rows','complete');
                ACC.Fc_may_ensmean(x,y,:) = tmp_corr(sub2ind(size(tmp_corr),[1:9],[1:9]));

        end
end

% Spatial masks
isebs = Grd.surveystrata_updated < 100 & isnan(squeeze(Hc.var(:,:,1)))==0;


%end

% Plot Figures

[X,Y] = mask2poly(Grd.lon_psi,Grd.lat_psi,isebs(2:end-1,2:end-1));

cmap = customcolormap(linspace(0,1,11), {'#523107','#523107','#bf812c','#e2c17e','#f3e9c4','#f6f4f4','#cae9e3','#81cdc1','#379692','#01665e','#003d2e'}); % brown-white-pool

figure(1); set(gcf, 'units','centimeters','position',[10 20 37 19]); set(gcf,'Color',[1 1 1]);
for m = 1:9
        ax(m) = subplot(3,3,m);
        m_proj('Mercator','long',[177 205],'lat',[53 66])
        fi = find(Grd.h > 1500); tmp_plot = squeeze(ACC.Fc_may_ensmean(:,:,m)); tmp_plot(fi) = nan;
        fi = find(tmp_plot < 0); tmp_plot(fi) = 0; % Mask points < 0 as 0 so that contourf doesn't mark them as white when caxis stops at 0
        [C,h] = m_contourf(Grd.lon_rho,Grd.lat_rho,tmp_plot,[0:.025:1]);
        set(h,'LineColor','none')
        hold on
        m_line(X,Y,'color','k','Linewidth',2)
        [C200,h200] = m_contour(Grd.lon_rho, Grd.lat_rho, Grd.h, [50 100 200], 'LineWidth', 0.5, 'Color', [0 0 0]);
        caxis([0 1])
        m_gshhs_i('patch',[0.5 0.5 0.5]);
        m_grid('xtick',[160:20:210],'ytick',[50:8:66],'linestyle','none','color','k','Fontsize',6);
        title([num2str(m),' month lead'],'Fontsize',12)
end
colormap(flipud(cmap))
cb = colorbar;
set(cb,'position',[ax(8).Position(1) ax(8).Position(2)-.052 ax(8).Position(3) .035],'orientation','horizontal')
xlabel(cb,'ACC','Fontsize',12)
set(gca,'SortMethod','ChildOrder') % To avoid warning message with export_fig





