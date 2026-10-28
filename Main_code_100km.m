clear;
close all;
folder = fileparts(matlab.desktop.editor.getActiveFilename); 
addpath(genpath(folder));

load 'conflict_data_100km';
load 'Conflict100km_n18316data_all.mat';

%% Regression modelling
%% Logistic regression variables processing
[Conflict100km_n18316data_nonan,shouldremove] = rmmissing(Conflict100km_n18316data_all);
% writetable(Conflict100km_n18316data_nonan,'Conflict100km_n12707.csv')

%%
num_conflict = sum(Conflict100km_n18316data_nonan.binary_conflict)

num_conflict_type1 = sum(Conflict100km_n18316data_nonan.conflict_type1>0)
num_conflict_type2 = sum(Conflict100km_n18316data_nonan.conflict_type2>0)
num_conflict_type3 = sum(Conflict100km_n18316data_nonan.conflict_type3>0)

%% niche or niche loss correlation with conflict density
conlfictdensity = R_GED24point2_9020_100km./pixelarea; 
conlfictdensity(isnan(conlfictdensity)) = 0;
log_conlfictdensity = log(conlfictdensity*1e4+1);

se_var = Conflict100km_n18316data_all.niche_comb_norm; % Conflict100km_n18316data_all.nicheloss_comb_lognorm; % or Conflict100km_n18316data_all.niche_comb_norm
max_se = max(se_var);
min_se = min(se_var);
segment = 7;
factorwindow = (max_se - min_se) /segment;
factorstep = factorwindow/2.0;
factorrange = [min_se max_se];
leastsample = 50;
runmean_vec  = runningmeandata(log_conlfictdensity(~shouldremove),  se_var(~shouldremove),  factorwindow, factorstep, factorrange, leastsample); 
binvec = factorrange(1):factorstep:factorrange(2)+factorstep/2.0;
plot(binvec, runmean_vec,'-','LineWidth', 2);  
axis square
set(gca, 'FontSize', 12);
set(gca, 'XDir', 'reverse');  
ylim([0 1])

% a = array2table([log_conlfictdensity(~shouldremove),Conflict100km_n18316data_all.niche_comb_norm(~shouldremove),Conflict100km_n18316data_all.nicheloss_comb_lognorm(~shouldremove)]);
% a.Properties.VariableNames = {'ln(Conflict density)','niche','niche loss'};
% writetable(a, 'E:\DATA\Conflict\Manuscript\final version\source data\figure source data.xlsx', 'Sheet', 'Figure 3', 'WriteMode', 'overwritesheet');

deleteExcelSheet('E:/DATA/Conflict/Manuscript/final version/source data/figure source data.xlsx', 'Figure2_part1');
deleteExcelSheet('E:/DATA/Conflict/Manuscript/final version/source data/figure source data.xlsx', 'Figure2_part2');
deleteExcelSheet('E:/DATA/Conflict/Manuscript/final version/source data/figure source data.xlsx', 'Figure2_part3');

%% Spearsman correlation test between socioeconomic variables and conflict density
HDI_log = log(HDI_1990_2015_100km + 1);                     % human development index  
GNI_log = log(gni_perCapita_1990_2021_100km + 1);           % gross national income per capita
gini = gini_disp_1990_2021_100km;                           % income GINI
GDPpercap_log = log(GDP_per_capita_1990_2020_100km_zonal_mean + 1); % GDP per capita 
langdiv_1dlog  = log(langdiv_1d_100km+1);

[rho, pval] = corr(HDI_log(~shouldremove), log_conlfictdensity(~shouldremove), 'Type', 'Spearman');
rho
pval
[rho, pval] = corr(GNI_log(~shouldremove), log_conlfictdensity(~shouldremove), 'Type', 'Spearman');
rho
pval
[rho, pval] = corr(gini(~shouldremove), log_conlfictdensity(~shouldremove), 'Type', 'Spearman');
rho
pval
[rho, pval] = corr(GDPpercap_log(~shouldremove), log_conlfictdensity(~shouldremove), 'Type', 'Spearman');
rho
pval
[rho, pval] = corr(langdiv_1dlog(~shouldremove), log_conlfictdensity(~shouldremove), 'Type', 'Spearman');
rho
pval
[rho, pval] = corr(Conflict100km_n18316data_all.nicheloss_comb_lognorm(~shouldremove), log_conlfictdensity(~shouldremove), 'Type', 'Spearman');
rho
pval