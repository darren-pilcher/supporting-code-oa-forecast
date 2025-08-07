% Script used to calculate the ACC scores for the seasonal forecasts
% Created by DJP 09/24

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

% Grab points during selected timeframe
ind_tmp = isbetween(var_t,datetime(years(1),01,01),datetime(years(end),12,31,23,59,59));

Hc.var = var_tmp(:,:,ind_tmp);
Hc.t = var_t(ind_tmp);
clear var*

% Eastern Bering Sea Bristol Bay mask and indices
BB = shaperead('./BBRKC_Shape_Files/BristolBay.shp');

% Convert Longitude to positive degrees East
BB.BoundingBox(:,1) = BB.BoundingBox(:,1) + 360;

BB.X = BB.X + 360;

% Make masks
[in_BB,on_BB] = inpolygon(Grd.lon_rho,Grd.lat_rho,BB.X,BB.Y);

isebs = Grd.surveystrata_updated < 100 & isnan(squeeze(Hc.var(:,:,1)))==0;
%isebs = (in_BB == 1 & isnan(squeeze(Hc.var(:,:,1)))==0); % This is BB mask, just uncomment here to get BB specific results

% Calculate EBS spatial average here
Hc.ebs = local(Hc.var,isebs,'weight',Grd.area_feast,'omitnan');

% Calculate monthly averages in prep for montly climatology
end_day = [31 28 31 30 31 30 31 31 30 31 30 31];

% Loop through each month
for m = 1:12
        day_st = datetime(years,m,1); day_en = datetime(years,m,end_day(m),23,59,59);
        for t = 1:length(years)
                tmp_ind = isbetween(Hc.t,day_st(t),day_en(t));
                Hc.ebs_mon(m,t) = mean(Hc.ebs(tmp_ind));
        end
end

% Load forecasts
% Load bottom pH and arag output from forecasts (Incorporate 1982-1983 later due to missing ensemble members)
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
					tmp_ebs = local(var_tmp,isebs,'weight',Grd.area_feast,'omitnan');
					Fc.ebs_mon_apr(e,m,t+2) = mean(tmp_ebs(tmp_ind));
				end	
			elseif i==2
				for m = 1:8
					day_st = datetime(years(t+2),m+4,1); day_en = datetime(years(t+2),m+4,end_day(m+4),23,59,59);
                                        tmp_ind = isbetween(tmp_t,day_st,day_en);
					tmp_ebs = local(var_tmp,isebs,'weight',Grd.area_feast,'omitnan');
                                        Fc.ebs_mon_may(e,m,t+2) = mean(tmp_ebs(tmp_ind));
				end	
				% a bit clunky but accounts for Jan of following year
				years_end = [1982:2011];
				day_st = datetime(years_end(t+3),1,1); day_en = datetime(years_end(t+3),1,end_day(1),23,59,59);
                                tmp_ind = isbetween(tmp_t,day_st,day_en);
                                tmp_ebs = local(var_tmp,isebs,'weight',Grd.area_feast,'omitnan');
                                Fc.ebs_mon_may(e,9,t+2) = mean(tmp_ebs(tmp_ind));	
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
                        tmp_ebs = local(var_tmp,isebs,'weight',Grd.area_feast,'omitnan');
			Fc.ebs_mon_may(e,m,t) = mean(tmp_ebs(tmp_ind));	
		end
		day_st = datetime(years(t+1),1,1); day_en = datetime(years(t+1),1,end_day(1),23,59,59);
                tmp_ind = isbetween(tmp_t,day_st,day_en);
                tmp_ebs = local(var_tmp,isebs,'weight',Grd.area_feast,'omitnan');
                Fc.ebs_mon_may(e,9,t) = mean(tmp_ebs(tmp_ind));
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
		tmp_ebs = local(var_tmp,isebs,'weight',Grd.area_feast,'omitnan');
                Fc.ebs_mon_apr(e,m,t) = mean(tmp_ebs(tmp_ind));
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
		tmp_ebs = local(var_tmp,isebs,'weight',Grd.area_feast,'omitnan');
                Fc.ebs_mon_apr(e,m,t) = mean(tmp_ebs(tmp_ind));
        end
end


% Fill in nans for missing ensemble members
Fc.ebs_mon_apr(3,:,1) = nan; Fc.ebs_mon_apr(2,:,2) = nan;

% Detrend hindcast.  Start by calculating trend over annual average values
if (detrend)
	tmp = squeeze(nanmean(Hc.ebs_mon,1));
	Hc.lin_trend = trend(tmp); clear tmp

% Now remove annual lin trend from monthly hindcast output
	for t = 1:length(years)
		for m = 1:12
			Hc.ebs_mon(m,t) = Hc.ebs_mon(m,t) - (Hc.lin_trend*t);
		end
	end

% Now detrend forecasts
	for t = 1:length(years)
		for e = 1:3
			for m = 1:9
				Fc.ebs_mon_apr(e,m,t) = Fc.ebs_mon_apr(e,m,t) - (Hc.lin_trend*t);
			end

			for m = 1:9
				Fc.ebs_mon_may(e,m,t) = Fc.ebs_mon_may(e,m,t) - (Hc.lin_trend*t);
			end
		end
	end
			

end % detrend

% Compute hindcast monthly anomalies
for t = 1:length(years)
	for m = 1:12
		Hc.ebs_mon_anom(m,t) = squeeze(Hc.ebs_mon(m,t)) - squeeze(nanmean(Hc.ebs_mon(m,:),2));
	end
end

% Compute forecast monthly anomalies
% April
for t = 1:length(years)
	for e = 1:3
		for m = 1:9
			Fc.ebs_mon_apr_anom(e,m,t) = squeeze(Fc.ebs_mon_apr(e,m,t)) - squeeze(nanmean(Fc.ebs_mon_apr(e,m,:),3));
		end
	end
% Also calculate ensemble mean anomalies
	for m = 1:9
		Fc.ebs_mon_apr_ensmean_anom(m,t) = squeeze(nanmean(Fc.ebs_mon_apr(:,m,t),1)) - squeeze(nanmean(nanmean(Fc.ebs_mon_apr(:,m,:),1),3));; 
	end

end

% May
for t = 1:length(years)
        for e = 1:3
                for m = 1:9
                        Fc.ebs_mon_may_anom(e,m,t) = squeeze(Fc.ebs_mon_may(e,m,t)) - squeeze(nanmean(Fc.ebs_mon_may(e,m,:),3));
                end
        end
% Also calculate ensemble mean anomalies
        for m = 1:9
                Fc.ebs_mon_may_ensmean_anom(m,t) = squeeze(nanmean(Fc.ebs_mon_may(:,m,t),1)) - squeeze(nanmean(nanmean(Fc.ebs_mon_may(:,m,:),1),3));;
        end
end

% Compute ACC between forecasts and hindcast
for e = 1:3
	tmp_hc = squeeze(Hc.ebs_mon_anom(4:12,:))'; tmp_fc = squeeze(Fc.ebs_mon_apr_anom(e,:,:))';
	tmp_corr = corr(tmp_hc,tmp_fc,'rows','complete');
	ACC.Fc_apr(e,:) = tmp_corr(sub2ind(size(tmp_corr),[1:9],[1:9])); 		

        tmp_hc = squeeze(Hc.ebs_mon_anom(5:12,:))'; tmp_fc = squeeze(Fc.ebs_mon_may_anom(e,:,:))';
	tmp_hc_end = [squeeze(Hc.ebs_mon_anom(1,2:end))';nan];	
	tmp_hc = [tmp_hc tmp_hc_end];
        tmp_corr = corr(tmp_hc,tmp_fc,'rows','complete');
        ACC.Fc_may(e,:) = tmp_corr(sub2ind(size(tmp_corr),[1:9],[1:9])); 

end		


% Same as above but for ensemble mean forecast
tmp_hc = squeeze(Hc.ebs_mon_anom(4:12,:))'; tmp_fc = Fc.ebs_mon_apr_ensmean_anom';
tmp_corr = corr(tmp_hc,tmp_fc,'rows','complete');
ACC.Fc_apr_ensmean = tmp_corr(sub2ind(size(tmp_corr),[1:9],[1:9]));

tmp_hc = squeeze(Hc.ebs_mon_anom(5:12,:))'; tmp_fc = Fc.ebs_mon_may_ensmean_anom';
tmp_hc_end = [squeeze(Hc.ebs_mon_anom(1,2:end))';nan];
tmp_hc = [tmp_hc tmp_hc_end];
tmp_corr = corr(tmp_hc,tmp_fc,'rows','complete');
ACC.Fc_may_ensmean = tmp_corr(sub2ind(size(tmp_corr),[1:9],[1:9]));

% Also calculate persistence forecast
tmp_hc = squeeze(Hc.ebs_mon_anom(4:12,:))'; 
tmp = squeeze(Hc.ebs_mon_anom(3,:))'; tmp_pc = repmat(tmp,[1 9]);
Pc.apr = tmp_pc' + Hc.ebs_mon(4:12,:); 
tmp_corr = corr(tmp_hc,tmp_pc,'rows','complete');
ACC.Pc_apr = tmp_corr(sub2ind(size(tmp_corr),[1:9],[1:9])); 

tmp_hc = squeeze(Hc.ebs_mon_anom(5:12,:))';
tmp_hc_end = [squeeze(Hc.ebs_mon_anom(1,2:end))';nan]; tmp_hc = [tmp_hc tmp_hc_end];
tmp = squeeze(Hc.ebs_mon_anom(4,:))'; tmp_pc = repmat(tmp,[1 9]);
Hc.ebs_mon_may = [Hc.ebs_mon(5:12,:);[Hc.ebs_mon(1,2:end)';nan]'];
Pc.may = tmp_pc' + Hc.ebs_mon_may;
tmp_corr = corr(tmp_hc,tmp_pc,'rows','complete');
ACC.Pc_may = tmp_corr(sub2ind(size(tmp_corr),[1:9],[1:9]));

Fc.ebs_mon_apr_ensmean = squeeze(nanmean(Fc.ebs_mon_apr,1));
Fc.ebs_mon_may_ensmean = squeeze(nanmean(Fc.ebs_mon_may,1));

% Calculate effective degrees of freedom following Bretherton et al., (1999) equation 30
% A lot of this code from  Kearney et al., 2021
% calculate autocorrelation across lags from 0 to 27 years

% April
na = length(years); nm = length(years)-1; % April vs May forecasts.  Cutoff last May year forecast because it doesn't have Jan of 2011
for m = 1:9
	denom = 0; % variable stores denominator of eq. 30
	for tau = 0:27        
		a = corrcoef(Hc.ebs_mon(m+3,1:(na-tau)),Hc.ebs_mon(m+3,(1+tau):na));
        	rho_X1 = a(1,2);
        	a = corrcoef(Fc.ebs_mon_apr_ensmean(m,1:(na-tau)),Fc.ebs_mon_apr_ensmean(m,(1+tau):na));
        	rho_Y1 = a(1,2);
        	denom = denom + (1 - tau/na)*rho_X1*rho_Y1;
	end
	Fc.n_eff_dyn_apr(m) = na/denom;

	% Now for persistence forecast
	denom = 0; % variable stores denominator of eq. 30
        for tau = 0:27
                a = corrcoef(Hc.ebs_mon(m+3,1:(na-tau)),Hc.ebs_mon(m+3,(1+tau):na));
                rho_X1 = a(1,2);
                a = corrcoef(Pc.apr(m,1:(na-tau)),Pc.apr(m,(1+tau):na));
                rho_Y1 = a(1,2);
                denom = denom + (1 - tau/na)*rho_X1*rho_Y1;
        end
        Pc.n_eff_dyn_apr(m) = na/denom;
% Now for May
        denom = 0; % variable stores denominator of eq. 30
        for tau = 0:26
                a = corrcoef(Hc.ebs_mon_may(m,1:(nm-tau)),Hc.ebs_mon_may(m,(1+tau):nm));
                rho_X1 = a(1,2);
                a = corrcoef(Fc.ebs_mon_may_ensmean(m,1:(nm-tau)),Fc.ebs_mon_may_ensmean(m,(1+tau):nm));
                rho_Y1 = a(1,2);
                denom = denom + (1 - tau/nm)*rho_X1*rho_Y1;
        end
        Fc.n_eff_dyn_may(m) = nm/denom;

        % Now for persistence forecast
        denom = 0; % variable stores denominator of eq. 30
        for tau = 0:26
                a = corrcoef(Hc.ebs_mon_may(m,1:(nm-tau)),Hc.ebs_mon_may(m,(1+tau):nm));
                rho_X1 = a(1,2);
                a = corrcoef(Pc.may(m,1:(nm-tau)),Pc.may(m,(1+tau):nm));
                rho_Y1 = a(1,2);
                denom = denom + (1 - tau/nm)*rho_X1*rho_Y1;
        end
        Pc.n_eff_dyn_may(m) = nm/denom;

end

% Test if ACC is signficantly greater than 0
for m = 1:9
	% transform to Z-statistic (normally distributed)
        Z_mean = 0.5*log( (1+ACC.Fc_apr_ensmean(m))/(1-ACC.Fc_apr_ensmean(m)) );
        % standard deviation of Z-statistic
        Z_sigma = 1/sqrt(Fc.n_eff_dyn_apr(m) - 3);
        % Calculate the 5% lower bound
        Z_alpha = norminv(0.05,Z_mean,Z_sigma);
        % Transform back into correlation space using inverse
        % Fisher Z (e.g., eq. 8 of Vecchi et al., 2013)
        A.r_low = (exp(2*Z_alpha) - 1)/(exp(2*Z_alpha) + 1);
        % store cases where correlation is significantly greater
        % than 0

	Fc.gt0_apr(m) = A.r_low > 0;


	% transform to Z-statistic (normally distributed)
        Z_mean = 0.5*log( (1+ACC.Fc_may_ensmean(m))/(1-ACC.Fc_may_ensmean(m)) );
        % standard deviation of Z-statistic
        Z_sigma = 1/sqrt(Fc.n_eff_dyn_may(m) - 3);
        % Calculate the 5% lower bound
        Z_alpha = norminv(0.05,Z_mean,Z_sigma);
        % Transform back into correlation space using inverse
        % Fisher Z (e.g., eq. 8 of Vecchi et al., 2013)
        A.r_low = (exp(2*Z_alpha) - 1)/(exp(2*Z_alpha) + 1);
        % store cases where correlation is significantly greater
        % than 0

        Fc.gt0_may(m) = A.r_low > 0;

% Also calculate for persistence

        % transform to Z-statistic (normally distributed)
        Z_mean = 0.5*log( (1+ACC.Pc_apr(m))/(1-ACC.Pc_apr(m)) );
        % standard deviation of Z-statistic
        Z_sigma = 1/sqrt(Pc.n_eff_dyn_apr(m) - 3);
        % Calculate the 5% lower bound
        Z_alpha = norminv(0.05,Z_mean,Z_sigma);
        % Transform back into correlation space using inverse
        % Fisher Z (e.g., eq. 8 of Vecchi et al., 2013)
        A.r_low = (exp(2*Z_alpha) - 1)/(exp(2*Z_alpha) + 1);
        % store cases where correlation is significantly greater
        % than 0

        Pc.gt0_apr(m) = A.r_low > 0;


        % transform to Z-statistic (normally distributed)
        Z_mean = 0.5*log( (1+ACC.Pc_may(m))/(1-ACC.Pc_may(m)) );
        % standard deviation of Z-statistic
        Z_sigma = 1/sqrt(Pc.n_eff_dyn_may(m) - 3);
        % Calculate the 5% lower bound
        Z_alpha = norminv(0.05,Z_mean,Z_sigma);
        % Transform back into correlation space using inverse
        % Fisher Z (e.g., eq. 8 of Vecchi et al., 2013)
        A.r_low = (exp(2*Z_alpha) - 1)/(exp(2*Z_alpha) + 1);
        % store cases where correlation is significantly greater
        % than 0

        Pc.gt0_may(m) = A.r_low > 0;


end

% betterprc = percentile uncertainty at which S1 is better than S2
% i.e. 5 -> better at 5% threshold, 10 -> better at 10% threshold
for m = 1:9
    % generate 1000 samples of dynamical forecast correlation
    Z_mean_dyn = 0.5*log( (1+ACC.Fc_apr_ensmean(m))/(1-ACC.Fc_apr_ensmean(m)) );
    Z_sigma_dyn = 1/sqrt(Fc.n_eff_dyn_apr(m)-3);
    Z_rand_dyn = normrnd(Z_mean_dyn*ones(1000,1),Z_sigma_dyn);

    % generate 1000 samples of persistence forecast correlation
    Z_mean_pers = 0.5*log( (1+ACC.Pc_apr(m))/(1-ACC.Pc_apr(m)) );
    Z_sigma_pers = 1/sqrt(Pc.n_eff_dyn_apr(m)-3);
    Z_rand_pers = normrnd(Z_mean_pers*ones(1000,1),Z_sigma_pers);

    % assign significance to dynamical-persistence
    aa = find(Z_rand_pers - Z_rand_dyn > 0);

    Fc.betterprc_apr(m) = (size(aa,1)/1000)*100;

    % generate 1000 samples of dynamical forecast correlation
    Z_mean_dyn = 0.5*log( (1+ACC.Fc_may_ensmean(m))/(1-ACC.Fc_may_ensmean(m)) );
    Z_sigma_dyn = 1/sqrt(Fc.n_eff_dyn_may(m)-3);
    Z_rand_dyn = normrnd(Z_mean_dyn*ones(1000,1),Z_sigma_dyn);

    % generate 1000 samples of persistence forecast correlation
    Z_mean_pers = 0.5*log( (1+ACC.Pc_may(m))/(1-ACC.Pc_may(m)) );
    Z_sigma_pers = 1/sqrt(Pc.n_eff_dyn_may(m)-3);
    Z_rand_pers = normrnd(Z_mean_pers*ones(1000,1),Z_sigma_pers);

    % assign significance to dynamical-persistence
    aa = find(Z_rand_pers - Z_rand_dyn > 0);

    Fc.betterprc_may(m) = (size(aa,1)/1000)*100;
end

% Use a 5% thresholds for dynamic over persistence
Fc.pc5_apr = Fc.betterprc_apr < 5;
Fc.pc5_may = Fc.betterprc_may < 5;

% Save .mat file with stats for generating grid in plot_forecast_skill_scores_grid.m script
arag_surf_stats.ACC(1,:) = ACC.Fc_apr_ensmean;
arag_surf_stats.ACC(2,:) = ACC.Fc_may_ensmean;
arag_surf_stats.gt0(1,:) = Fc.gt0_apr;
arag_surf_stats.gt0(2,:) = Fc.gt0_may;
arag_surf_stats.pc5(1,:) = Fc.pc5_apr;
arag_surf_stats.pc5(2,:) = Fc.pc5_may;

Pc_pH_surf_stats.ACC(1,:) = ACC.Pc_apr;
Pc_pH_surf_stats.ACC(2,:) = ACC.Pc_may;
Pc_pH_surf_stats.gt0(1,:) = Pc.gt0_apr;
Pc_pH_surf_stats.gt0(2,:) = Pc.gt0_may;

%save('./arag_surf_stats.mat','arag_surf_stats','-v7.3')
%save('./Pc_pH_surf_stats.mat','Pc_pH_surf_stats','-v7.3')

% Plot Figures
figure(1); set(gcf, 'units','centimeters','position',[10 40 19 15]); set(gcf,'Color',[1 1 1]);
subplot(2,1,1)
plot(ACC.Fc_apr','Linewidth',2,'Marker','o')
hold on
plot(mean(ACC.Fc_apr,1),'k-o','Linewidth',2)
plot(ACC.Fc_apr_ensmean,'Color',[192/255 192/255 192/255],'Marker','o','Linewidth',2)
%plot(mean(ACC.Fc_apr_ebs,1),'k-o','Linewidth',2)
plot(ones(1,9)*0.5,'k--','Linewidth',1)
%plot(ones(1,9)*0.6,'k--','Linewidth',1)
axis([1 9 0 1])
ylabel('ACC','Fontsize',14)
%title(['April Init ',VAR,' ',depth,' Forecasts'],'Fontsize',16)
title('April Initialized pH Bot Forecasts','Fontsize',16)
lgd = legend('Ens 1','Ens 2','Ens 3','Mean of Ens','Ens mean','Location','Southwest');
lgd.Location = "south"; lgd.Orientation = "horizontal";

subplot(2,1,2)
plot(ACC.Fc_may','Linewidth',2,'Marker','o')
hold on
plot(mean(ACC.Fc_may,1),'k-o','Linewidth',2)
plot(ACC.Fc_may_ensmean,'Color',[192/255 192/255 192/255],'Marker','o','Linewidth',2)
%plot(mean(ACC.Fc_may_ebs,1),'k-o','Linewidth',2)
plot(ones(1,9)*0.5,'k--','Linewidth',1)
%plot(ones(1,9)*0.6,'k--','Linewidth',1)
axis([1 9 0 1])
xlabel('Lead time (months)','Fontsize',14)
ylabel('ACC','Fontsize',14)
%title(['May Init ',VAR,' ',depth,' Forecasts'],'Fontsize',16)
title('May Initialized pH Bot Forecasts','Fontsize',16)
lgd = legend('Ens 1','Ens 2','Ens 3','Mean of Ens','Ens mean','Location','Southwest');
lgd.Location = "south"; lgd.Orientation = "horizontal";

%export_fig('/gscratch/cicoes/GR011909_oap/pilchd/Bering10k/figures/forecasts_jun2022/pH_bot_acctimeseries_allmembers','-png','tif','jpg')
%f = gcf;
%exportgraphics(f,'/gscratch/cicoes/GR011909_oap/pilchd/Bering10k/figures/forecasts_jun2022/test.pdf');

figure(2); set(gcf, 'units','centimeters','position',[10 45 19 15]); set(gcf,'Color',[1 1 1]);
subplot(2,1,1)
plot(ACC.Fc_apr_ensmean,'k','Linewidth',2)
hold on
plot(ACC.Pc_apr,'Color',[192/255 192/255 192/255],'Linewidth',2)
plot(ones(1,9)*0.5,'k--','Linewidth',1)
tmp = ACC.Fc_apr_ensmean; tmp(~Fc.gt0_apr)=nan;
plot(tmp,'k','Marker','o','MarkerFaceColor','k','LineStyle','none')
tmp = ACC.Fc_apr_ensmean; tmp(Fc.gt0_apr)=nan;
plot(tmp,'k','Marker','o','LineStyle','none')
tmp = ACC.Pc_apr; tmp(~Pc.gt0_apr)=nan;
plot(tmp,'Color',[192/255 192/255 192/255],'Marker','o','MarkerFaceColor',[192/255 192/255 192/255],'LineStyle','none')
tmp = ACC.Pc_apr; tmp(Pc.gt0_apr)=nan;
plot(tmp,'Color',[192/255 192/255 192/255],'Marker','o','LineStyle','none')
% Plot forecast > persistence markers above forecast line
tmp = ACC.Fc_apr_ensmean - 0.04; tmp(~Fc.pc5_apr) = nan;
plot(tmp,'k','Marker','*','LineStyle','none')
axis([1 9 -.05 1])
ylabel('ACC','Fontsize',14)
%title(['April Init ',VAR,' ',depth,' Forecasts'],'Fontsize',16)
title('April Initialized Bottom Water pH Forecasts','Fontsize',16)
legend('Dynamic','Persistence','Location','Southeast')

subplot(2,1,2)
plot(ACC.Fc_may_ensmean,'k','Linewidth',2)
hold on
plot(ACC.Pc_may,'Color',[192/255 192/255 192/255],'Linewidth',2)
plot(ones(1,9)*0.5,'k--','Linewidth',1)
tmp = ACC.Fc_may_ensmean; tmp(~Fc.gt0_may)=nan;
plot(tmp,'k','Marker','o','MarkerFaceColor','k','LineStyle','none')
tmp = ACC.Fc_may_ensmean; tmp(Fc.gt0_may)=nan;
plot(tmp,'k','Marker','o','LineStyle','none')
tmp = ACC.Pc_may; tmp(~Pc.gt0_may)=nan;
plot(tmp,'Color',[192/255 192/255 192/255],'Marker','o','MarkerFaceColor',[192/255 192/255 192/255],'LineStyle','none')
tmp = ACC.Pc_may; tmp(Pc.gt0_may)=nan;
plot(tmp,'Color',[192/255 192/255 192/255],'Marker','o','LineStyle','none')
tmp = ACC.Fc_may_ensmean - 0.04; tmp(~Fc.pc5_may) = nan;
plot(tmp,'k','Marker','*','LineStyle','none')
axis([1 9 -.05 1])
ylabel('ACC','Fontsize',14)
%title(['May Init ',VAR,' ',depth,' Forecasts'],'Fontsize',16)
title('May Initialized Bottom Water pH Forecasts','Fontsize',16)
legend('Dynamic','Persistence','Location','Southeast')

%export_fig('/gscratch/cicoes/GR011909_oap/pilchd/Bering10k/figures/forecasts_jun2022/pH_surface_acc_ensmean_timeseries_sig','-pdf')
%f = gcf;
%exportgraphics(f,'/gscratch/cicoes/GR011909_oap/pilchd/Bering10k/figures/forecasts_jun2022/arag_bot_acc_ensmean_timeseries_sig_detrend.pdf')

figure(3); set(gcf, 'units','centimeters','position',[10 10 19 11]); set(gcf,'Color',[1 1 1]);
[p1] = plot(years,squeeze(nanmean(Hc.ebs_mon,1)),'k-o','Linewidth',2);
hold on
[p2] = polyplot(years,squeeze(nanmean(Hc.ebs_mon,1)),'k--','Linewidth',2);
set(gca,'xtick',[1982:4:2010],'Fontsize',12)
%set(gca,'ytick',[1.4:.1:1.9],'Fontsize',12)
%ylabel(VAR,'Fontsize',16)
axis([1982 2010 7.88 8])
ylabel('{\Omega_{arag}}','Fontsize',16)
xlabel('Year','Fontsize',16)
%legend([p4,p5,p6],['Slope = ',num2str(trend(var_shelf_mean_full),2)])
%title(['Bering Sea Shelf ',t_label,' ',depth,' ',VAR],'Fontsize',16)
legend([p2],['Trend = ',num2str(trend(squeeze(nanmean(Hc.ebs_mon,1))),2)],'Location','Northeast')
title('Bering Sea Shelf Annual Bottom Water pH','Fontsize',18)

%export_fig('/gscratch/jisao/pilchd/Bering10k/figures/forecasts_jun2022/pH_surf_hindcast_trend','-pdf')
%f = gcf;
%exportgraphics(f,'/gscratch/cicoes/GR011909_oap/pilchd/Bering10k/figures/forecasts_jun2022/arag_bot_hindcast_trend.pdf')

