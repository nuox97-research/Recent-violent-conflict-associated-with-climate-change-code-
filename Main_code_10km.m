clear;
folder = fileparts(matlab.desktop.editor.getActiveFilename); 
addpath(genpath(folder));

load 'conflict_data_10km';

%% temperature niche 10 km （main figure 1ac）

popd_2005AD_10km(isnan(popd_2005AD_10km)) = 0;
popd_1980AD_10km(isnan(popd_1980AD_10km)) = 0;

mat_range   = [-5.0 35.0]; 
leastsample = 100;
C = linspecer(9,'sequential'); % line color setting function

popd_min = 0; % minimum population density
fltr = ~isnan(ERA5_MAT_1960_1990_10km) & ~isnan(GLDAS_MAT_1960_1990_10km) & ~isnan(CRU_MAT_1960_1990_10km) & popd_2005AD_10km > popd_min; % filter out population density = 0 (only applied to running mean plot of confict density)
mat6090_popfltr  = [ERA5_MAT_1960_1990_10km(fltr) CRU_MAT_1960_1990_10km(fltr) GLDAS_MAT_1960_1990_10km(fltr) ];
mat9020_popfltr  = [ERA5_MAT_1990_2020_10km(fltr) CRU_MAT_1990_2020_10km(fltr) GLDAS_MAT_1990_2020_10km(fltr) ];
fltr2 = ~isnan(ERA5_MAT_1960_1990_10km) & ~isnan(GLDAS_MAT_1960_1990_10km) & ~isnan(CRU_MAT_1960_1990_10km);
mat9020_fltr2  = [ERA5_MAT_1990_2020_10km(fltr2) CRU_MAT_1990_2020_10km(fltr2) GLDAS_MAT_1990_2020_10km(fltr2) ];
mat6090_fltr2  = [ERA5_MAT_1960_1990_10km(fltr2) CRU_MAT_1960_1990_10km(fltr2) GLDAS_MAT_1960_1990_10km(fltr2) ];
smallbaseno = 0.005; % avoiding Inf (divided by zero)

% temperature niche and niche loss
nichechange_era_10km = n_t_1_era_10km - n_t_0_era_10km;
niche1dloss10km_era = nichechange_era_10km; 
niche1dloss10km_era(niche1dloss10km_era>0) = 0;
niche1dloss10km_era(isnan(niche1dloss10km_era)) = 0;
relativeloss_era = niche1dloss10km_era./(n_t_0_era_10km + smallbaseno);

nichechange_cru_10km = n_t_1_cru_10km - n_t_0_cru_10km;
niche1dloss10km_cru = nichechange_cru_10km; 
niche1dloss10km_cru(niche1dloss10km_cru>0) = 0;
niche1dloss10km_cru(isnan(niche1dloss10km_cru)) = 0;
relativeloss_cru = niche1dloss10km_cru./(n_t_0_cru_10km + smallbaseno);

nichechange_gld_10km = n_t_1_gld_10km - n_t_0_gld_10km;
niche1dloss10km_gld = nichechange_gld_10km; 
niche1dloss10km_gld(niche1dloss10km_gld>0) = 0;
niche1dloss10km_gld(isnan(niche1dloss10km_gld)) = 0;
relativeloss_gld = niche1dloss10km_gld./(n_t_0_gld_10km + smallbaseno);

% human & conflict temperature niche
conflictfreq10km = R_GED24point2_9020_10km./pixelarea; 
conflictfreq10km(isnan(conflictfreq10km)) = 0;

scalefactor = 1.2; % a scale factor to plot together conflict/human/niche (with different units) running mean curves   
HCN_mat_rm(popd_1980AD_10km(fltr2) , mat6090_fltr2,  mat_range , leastsample, 2, scalefactor, 1); 
scalefactor = 0.5;
HCN_mat_rm(conflictfreq10km(fltr) , mat9020_popfltr,  mat_range , leastsample, 9, scalefactor, 1); 

% relative niche loss
scalefactor = 13;
rm_vec_rel_era  = runningmeandata(relativeloss_era(fltr2),  ERA5_MAT_1990_2020_10km(fltr2),  2.0, 1.0, mat_range,leastsample);
rm_vec_rel_cru  = runningmeandata(relativeloss_cru(fltr2),  CRU_MAT_1990_2020_10km(fltr2),   2.0, 1.0, mat_range,leastsample);
rm_vec_rel_gld  = runningmeandata(relativeloss_gld(fltr2),  GLDAS_MAT_1990_2020_10km(fltr2), 2.0, 1.0, mat_range,leastsample);

rm_rel_vec = [rm_vec_rel_era; rm_vec_rel_cru; rm_vec_rel_gld];
maxcur_rel_out = prctile(rm_rel_vec,95)  * scalefactor;                  % 95 percentile of running mean ensemble 
mincur_rel_out = prctile(rm_rel_vec,5)   * scalefactor;                  %  5 percentile of running mean ensemble 
meancurRM_rel_out = mean(rm_rel_vec,'omitnan'); 

binvec = mat_range(1):1.0:mat_range(2)+0.5;
fltout = ~isnan(meancurRM_rel_out);
h = fill([binvec(fltout),fliplr(binvec(fltout))],[mincur_rel_out(fltout),fliplr(maxcur_rel_out(fltout))],C(7,:));
set(h,'LineStyle','none','facealpha',0.3);
hold on;
plot(binvec,meancurRM_rel_out * scalefactor,'-','LineWidth',2, 'Color',C(7,:));
xxvec = -5:0.1:35; yyvec = xxvec*0;
plot(xxvec,yyvec,'-k');
hold off;
xlim([0 30]); ylim([-2.2 2.2]);
set(gcf, 'Position', [0, 0, 500, 700]);
title('Temperature niche & niche loss 10 km','FontSize',10);

relativeloss_tmpniche = [relativeloss_era  relativeloss_cru  relativeloss_gld];
% absltloss_tmpniche    = [niche1dloss10km_era  niche1dloss10km_cru  niche1dloss10km_gld];

%% MAT niche curve fitting （extended data figure 1a）

popuse = popd_1980AD_10km(fltr2);
[~,niche1dfunc6090] = niche1D(mean(mat6090_fltr2,2,'omitnan'), popuse);
[~,niche1dfunc9020] = niche1D(mean(mat9020_fltr2,2,'omitnan'), popuse);

matvec = -5:0.1:35;
matniche6090 = niche1dfunc6090(matvec);
matniche9020 = niche1dfunc9020(matvec);

figure; 
hold on;
plot(matvec,matniche6090,'--','Color',C(2,:),'Linewidth',2);

scalefactor = 1;
HCN_mat_rm(popd_1980AD_10km(fltr2) , mat6090_fltr2,  mat_range , leastsample, 2, scalefactor, 1); 
HCN_mat_rm(popd_1980AD_10km(fltr2) , mat9020_fltr2,  mat_range , leastsample, 1, scalefactor, 1); 
ylabel('Population density',FontSize=16);
set(gcf, 'Position', [0, 0, 500, 400]);
hold off;


%% precipitation niche 10 km（main figure 1bd）

map_range   = [0 4000]; 
leastsample = 100;
C = linspecer(9,'sequential'); 

popd_min = 0;
fltr = ~isnan(ERA5_MAP_1990_2020_10km) & ~isnan(CRU_MAP_1990_2020_10km) & ~isnan(GLDAS_MAP_1990_2020_10km) & ~isnan(Terra_MAP_1990_2020_10km) & popd_2005AD_10km > popd_min; % filter out population density = 0 (only applied to running mean plot of confict density)
map6090_popfltr  = [ERA5_MAP_1960_1990_10km(fltr) CRU_MAP_1960_1990_10km(fltr) GLDAS_MAP_1960_1990_10km(fltr) Terra_MAP_1960_1990_10km(fltr)];  
map9020_popfltr  = [ERA5_MAP_1990_2020_10km(fltr) CRU_MAP_1990_2020_10km(fltr) GLDAS_MAP_1990_2020_10km(fltr) Terra_MAP_1990_2020_10km(fltr)];
fltr2 = ~isnan(ERA5_MAP_1990_2020_10km) & ~isnan(CRU_MAP_1990_2020_10km) & ~isnan(GLDAS_MAP_1990_2020_10km) & ~isnan(Terra_MAP_1990_2020_10km);
map9020_fltr2  = [ERA5_MAP_1990_2020_10km(fltr2) CRU_MAP_1990_2020_10km(fltr2) GLDAS_MAP_1990_2020_10km(fltr2) Terra_MAP_1990_2020_10km(fltr2)];
map6090_fltr2  = [ERA5_MAP_1960_1990_10km(fltr2) CRU_MAP_1960_1990_10km(fltr2) GLDAS_MAP_1960_1990_10km(fltr2) Terra_MAP_1960_1990_10km(fltr2)];  
smallbaseno = 0.005;  % avoiding Inf (divided by zero)

% precipitation niche and niche loss
nichechange_era_10km = n_p_1_era_10km - n_p_0_era_10km;
niche1dloss10km_era = nichechange_era_10km; 
niche1dloss10km_era(niche1dloss10km_era>0) = 0;
niche1dloss10km_era(isnan(niche1dloss10km_era)) = 0;
relativeloss_era = niche1dloss10km_era./(n_p_0_era_10km + smallbaseno);

nichechange_cru_10km = n_p_1_cru_10km - n_p_0_cru_10km;
niche1dloss10km_cru = nichechange_cru_10km; 
niche1dloss10km_cru(niche1dloss10km_cru>0) = 0;
niche1dloss10km_cru(isnan(niche1dloss10km_cru)) = 0;
relativeloss_cru = niche1dloss10km_cru./(n_p_0_cru_10km + smallbaseno);

nichechange_gld_10km = n_p_1_gld_10km - n_p_0_gld_10km;
niche1dloss10km_gld = nichechange_gld_10km; 
niche1dloss10km_gld(niche1dloss10km_gld>0) = 0;
niche1dloss10km_gld(isnan(niche1dloss10km_gld)) = 0;
relativeloss_gld = niche1dloss10km_gld./(n_p_0_gld_10km + smallbaseno);

nichechange_ter_10km = n_p_1_ter_10km - n_p_0_ter_10km;
niche1dloss10km_ter = nichechange_ter_10km; 
niche1dloss10km_ter(niche1dloss10km_ter>0) = 0;
niche1dloss10km_ter(isnan(niche1dloss10km_ter)) = 0;
relativeloss_ter = niche1dloss10km_ter./(n_p_0_gld_10km + smallbaseno);

% human & conflict precipitation niche
conflictfreq10km = R_GED24point2_9020_10km./pixelarea; 
conflictfreq10km(isnan(conflictfreq10km)) = 0;
scalefactor = 1.7; % a scale factor to plot together conflict/human/niche (with different units) running mean curves   
HCN_map_rm(popd_1980AD_10km(fltr2),  map6090_fltr2,  map_range, leastsample, 2, scalefactor, 1); 
scalefactor = 1.1;
HCN_map_rm(conflictfreq10km(fltr), map9020_popfltr,  map_range, leastsample, 9, scalefactor, 1); 

% relative niche loss
rm_vec_era_relat  = runningmeandata(relativeloss_era(fltr2),  ERA5_MAP_1990_2020_10km(fltr2),  200, 100, map_range,leastsample);
rm_vec_cru_relat  = runningmeandata(relativeloss_cru(fltr2),  CRU_MAP_1990_2020_10km(fltr2),   200, 100, map_range,leastsample);
rm_vec_gld_relat  = runningmeandata(relativeloss_gld(fltr2),  GLDAS_MAP_1990_2020_10km(fltr2), 200, 100, map_range,leastsample);
rm_vec_ter_relat  = runningmeandata(relativeloss_ter(fltr2),  Terra_MAP_1990_2020_10km(fltr2), 200, 100, map_range,leastsample);

scalefactor = 12; 
rm_vec_relat = [rm_vec_era_relat; rm_vec_cru_relat; rm_vec_gld_relat; rm_vec_ter_relat];
maxcur_out_relat = prctile(rm_vec_relat,95)  * scalefactor;                  % 95 percentile of running mean ensemble 
mincur_out_relat = prctile(rm_vec_relat,5)   * scalefactor;                  %  5 percentile of running mean ensemble 
meancurRM_out_relat = mean(rm_vec_relat,'omitnan'); 

binvec = map_range(1):100:map_range(2)+50;
fltout = ~isnan(meancurRM_out_relat);
h = fill([binvec(fltout),fliplr(binvec(fltout))],[mincur_out_relat(fltout),fliplr(maxcur_out_relat(fltout))],C(7,:));
set(h,'LineStyle','none','facealpha',0.3);
plot(binvec,meancurRM_out_relat * scalefactor,'-','LineWidth',2, 'Color',C(7,:));
xxvec = 0:1:4000; yyvec = xxvec*0;
plot(xxvec,yyvec,'-k');
hold off;
xlim([0 3500]); ylim([-2.2 2.2]);
set(gcf, 'Position', [200, 200, 500, 700]);
title('precipitation niche & niche loss 10 km','FontSize',10);

relativeloss_preniche = [relativeloss_era relativeloss_cru relativeloss_gld relativeloss_ter];
absltloss_preniche    = [niche1dloss10km_era  niche1dloss10km_cru  niche1dloss10km_gld  niche1dloss10km_ter];


%% MAP niche curve fitting (extended data figure 1b)

popuse = popd_1980AD_10km(fltr2);
[~,niche1dfunc6090] = niche1D_map(mean(map6090_fltr2,2,'omitnan'), popuse);
[~,niche1dfunc9020] = niche1D_map(mean(map9020_fltr2,2,'omitnan'), popuse);

mapvec = 0:10:4000;
mapniche6090 = niche1dfunc6090(log(mapvec+1));
mapniche9020 = niche1dfunc9020(log(mapvec+1));

figure; hold on;
plot(mapvec,mapniche6090,'--','Color',C(2,:),'Linewidth',2);
scalefactor = 1;
HCN_map_rm(popd_1980AD_10km(fltr2),  map6090_fltr2,  map_range, leastsample, 2, scalefactor, 1); 
HCN_map_rm(popd_1980AD_10km(fltr2),  map9020_fltr2,  map_range, leastsample, 1, scalefactor, 1); 
ylabel('Population density',FontSize=16);
set(gcf, 'Position', [0, 0, 500, 400]);


%% Human & conflict 2D niche 10 km (main figure 2)
%% Human 2D niche

map_range = [0 4000];                               % data range mean annual temperature
mat_range = [-10 40];                               % data range mean annual precipitation
datarange =  [map_range(1) map_range(2) mat_range(1) mat_range(2)];
out_res = 100; 
nichefactor = 1e6;

mat9020  = [ERA5_MAT_1990_2020_10km  CRU_MAT_1990_2020_10km  GLDAS_MAT_1990_2020_10km];
map9020  = [ERA5_MAP_1990_2020_10km  CRU_MAP_1990_2020_10km  GLDAS_MAP_1990_2020_10km  Terra_MAP_1990_2020_10km];
mat9020_use = mean(mat9020, 2, 'omitnan');
map9020_use = mean(map9020, 2, 'omitnan');

mat6090  = [ERA5_MAT_1960_1990_10km  CRU_MAT_1960_1990_10km  GLDAS_MAT_1960_1990_10km];
map6090  = [ERA5_MAP_1960_1990_10km  CRU_MAP_1960_1990_10km  GLDAS_MAP_1960_1990_10km  Terra_MAP_1960_1990_10km];
mat6090_use = mean(mat6090, 2, 'omitnan');
map6090_use = mean(map6090, 2, 'omitnan');

surf_smooth = 7;     % default smoothing parameter
heatmap_smooth = 7;  % default smoothing parameter

% 2D niche with smoothed surface
[Surf_niche_ref] = niche2D_matmap(mat6090_use, map6090_use, popd_1980AD_10km, 1, out_res, surf_smooth, heatmap_smooth); 
title('reference niche');    
ylim([-5 35]);
figure;

% contour of human climate niche
xx = linspace(map_range(1),map_range(2),out_res);                 
yy = linspace(mat_range(1),mat_range(2),out_res);
[outx, outy] = meshgrid(xx,yy);

[~,niche_contour] = contour(outx,outy,Surf_niche_ref');
maxnicheval = max(Surf_niche_ref,[],'all');
contour_level_vec = [maxnicheval*0.9  maxnicheval*0.7  maxnicheval*0.5  maxnicheval*0.3  maxnicheval*0.1];
niche_contour.LevelList = contour_level_vec;
colormap('jet'); axis square;

Niche_contourTable = getContourLineCoordinates(niche_contour);
nichelevelvec = Niche_contourTable.Level;
nichelevelgroup = Niche_contourTable.Group;
nichecontour_x = Niche_contourTable.X;
nichecontour_y = Niche_contourTable.Y;
nichecontour_z = Niche_contourTable.Y * 0 + 100; % shown at top layer

%% Confilct niche

conflictfreq10km = R_GED24point2_9020_10km./pixelarea; 
conflictfreq10km(isnan(conflictfreq10km)) = 0;

% condlict density per 10 km grid cells 
surf_smooth = 7;
heatmap_smooth = 7;
[Surf_niche_conflicts] = niche2D_matmap(mat9020_use, map9020_use, conflictfreq10km./pixelarea, 1, out_res, surf_smooth, heatmap_smooth); 
title('Conflict niche','FontSize', 10);    
ylim([-5 35]);


hold on;
contour_color = colormap('Copper');
plot3(nichecontour_x(nichelevelgroup == 10),nichecontour_y(nichelevelgroup == 10),nichecontour_z(nichelevelgroup == 10),'Color',contour_color(1,:)  ,'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 9), nichecontour_y(nichelevelgroup == 9), nichecontour_z(nichelevelgroup == 9), 'Color',contour_color(10,:) ,'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 8), nichecontour_y(nichelevelgroup == 8), nichecontour_z(nichelevelgroup == 8), 'Color',contour_color(20,:) ,'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 7), nichecontour_y(nichelevelgroup == 7), nichecontour_z(nichelevelgroup == 7), 'Color',contour_color(60,:) ,'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 6), nichecontour_y(nichelevelgroup == 6), nichecontour_z(nichelevelgroup == 6), 'Color',contour_color(100,:),'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 5), nichecontour_y(nichelevelgroup == 5), nichecontour_z(nichelevelgroup == 5), 'Color',contour_color(140,:),'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 4), nichecontour_y(nichelevelgroup == 4), nichecontour_z(nichelevelgroup == 4), 'Color',contour_color(180,:),'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 3), nichecontour_y(nichelevelgroup == 3), nichecontour_z(nichelevelgroup == 3), 'Color',contour_color(200,:),'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 2), nichecontour_y(nichelevelgroup == 2), nichecontour_z(nichelevelgroup == 2), 'Color',contour_color(220,:),'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 1), nichecontour_y(nichelevelgroup == 1), nichecontour_z(nichelevelgroup == 1), 'Color',contour_color(250,:),'LineWidth', 1.5);

[~,niche_contour] = contour(outx,outy,Surf_niche_conflicts');
maxnicheval = max(Surf_niche_conflicts,[],'all');
contour_level_vec = [maxnicheval*0.9  maxnicheval*0.8  maxnicheval*0.7  maxnicheval*0.6  maxnicheval*0.5  maxnicheval*0.4];
niche_contour.LevelList = contour_level_vec;

matcountry = [14.8  17.6  22.5];  % Afgan(n = 36730)  Syria(n = 83752)  
mapcountry = [357   340   290 ];
zcountry   = [100   100   100];

plot3(mapcountry(1:2), matcountry(1:2), zcountry(1:2),'ow','MarkerSize',10,'LineWidth',2 );

mycolor = colormap('turbo');    
colormap(mycolor(15:end,:));
ylim([-5 35]); xlim([0 4000]);
axis square;    xlabel("MAP (mm)"); ylabel ("MAT (^{o}C)"); title(''); 
set(gca, 'Fontname', 'Arial', 'Fontsize', 38);
xlabel('MAP (mm)', 'Fontsize', 41);
ylabel('MAT (^{o}C)', 'Fontsize', 41);
set(gcf,'Position',[1 49 1920 954]);


%% Vector plot overlain on the conflict niche heatmap 

% Focus  on areas where conflicts have occurred.
conflictfreq10km = R_GED24point2_9020_10km./pixelarea; 
conflictfreq10km(isnan(conflictfreq10km)) = 0;

[conflictfreq10km_pos_loc,~] = find(conflictfreq10km>0);
conflictfreq10km_pos = conflictfreq10km(conflictfreq10km_pos_loc);
mat9020_use_pos = mat9020_use(conflictfreq10km_pos_loc);
map9020_use_pos = map9020_use(conflictfreq10km_pos_loc);

% get points within the niche contour
num_bins = 50;  % num of bins
mat_edges = linspace(mat_range(1),mat_range(2), num_bins+1);
map_edges = linspace(map_range(1),map_range(2), num_bins+1);

mat_centers = (mat_edges(1:end-1) + mat_edges(2:end)) / 2;
map_centers = (map_edges(1:end-1) + map_edges(2:end)) / 2;
[Px, Ty] = meshgrid(map_centers,mat_centers);

[Surf_niche_conflicts_vectorrange] = niche2D_matmap(mat9020_use_pos, map9020_use_pos, conflictfreq10km_pos, 1, out_res, surf_smooth, heatmap_smooth); 
pden_conflict_con = Surf_niche_conflicts_vectorrange';
pden_conflict_con(pden_conflict_con>0) = 1;
pden_conflict_con(pden_conflict_con<=0) = -1;
[~, h] = contour(outx, outy, pden_conflict_con,[0 0],'LineStyle','none'); 
Niche_contourTable_vector = getContourLineCoordinates(h);
Niche_contourTable_vector = Niche_contourTable_vector(Niche_contourTable_vector.Group == 1,:);

inside_points = [];
x = Niche_contourTable_vector.X;
y = Niche_contourTable_vector.Y;
polygons = [x, y];

is_inside = false(length(Px(:)),1);
pxall = Px(:);
tyall = Ty(:);
for k = 1:length(Px(:))
    xq = pxall(k);
    yq = tyall(k);

    [in, ~] = inpolygon(xq, yq, polygons(:,1), polygons(:,2));
    if in
        is_inside(k,1) = true;
    end
end
is_inside_squ = reshape(is_inside,num_bins,num_bins);
inside_points = [inside_points; Px(is_inside_squ), Ty(is_inside_squ)];

% condlict density per 10 km grid cells 
surf_smooth = 7;
heatmap_smooth = 7;
[Surf_niche_conflicts] = niche2D_matmap(mat9020_use, map9020_use, conflictfreq10km, 1, out_res, surf_smooth, heatmap_smooth); 
title('Conflict niche','FontSize', 10);     
ylim([-5 35]);

hold on;
contour_color = colormap('Copper');
plot3(nichecontour_x(nichelevelgroup == 10),nichecontour_y(nichelevelgroup == 10),nichecontour_z(nichelevelgroup == 10),'Color',contour_color(1,:)  ,'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 9), nichecontour_y(nichelevelgroup == 9), nichecontour_z(nichelevelgroup == 9), 'Color',contour_color(10,:) ,'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 8), nichecontour_y(nichelevelgroup == 8), nichecontour_z(nichelevelgroup == 8), 'Color',contour_color(20,:) ,'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 7), nichecontour_y(nichelevelgroup == 7), nichecontour_z(nichelevelgroup == 7), 'Color',contour_color(60,:) ,'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 6), nichecontour_y(nichelevelgroup == 6), nichecontour_z(nichelevelgroup == 6), 'Color',contour_color(100,:),'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 5), nichecontour_y(nichelevelgroup == 5), nichecontour_z(nichelevelgroup == 5), 'Color',contour_color(140,:),'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 4), nichecontour_y(nichelevelgroup == 4), nichecontour_z(nichelevelgroup == 4), 'Color',contour_color(180,:),'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 3), nichecontour_y(nichelevelgroup == 3), nichecontour_z(nichelevelgroup == 3), 'Color',contour_color(200,:),'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 2), nichecontour_y(nichelevelgroup == 2), nichecontour_z(nichelevelgroup == 2), 'Color',contour_color(220,:),'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 1), nichecontour_y(nichelevelgroup == 1), nichecontour_z(nichelevelgroup == 1), 'Color',contour_color(250,:),'LineWidth', 1.5);

mycolor = colormap('turbo');    
colormap(mycolor(15:end,:))
ylim([-5 35]); 
xlim([0 4000]);
% axis square;    
xlabel("MAP (mm)"); ylabel ("MAT (^{o}C)");  title(''); 

% mat map data
idx = (isnan(mat6090_use) | isnan(mat9020_use) | isnan(map6090_use) | isnan(map9020_use) | ...
    mat6090_use<mat_range(1) | mat6090_use>mat_range(2) | mat9020_use<mat_range(1) | mat9020_use>mat_range(2) |...
    map6090_use<map_range(1) | map6090_use>map_range(2) | map9020_use<map_range(1) | map9020_use>map_range(2));
loc_idx = logical(zeros(size(idx)));
loc_idx(conflictfreq10km_pos_loc) = true;
idx = (idx & ~loc_idx);
mat6090_use_pos = mat6090_use(~idx);
map6090_use_pos = map6090_use(~idx);
mat9020_use_pos = mat9020_use(~idx);
map9020_use_pos = map9020_use(~idx);

% calculate the arrow direction
[~, ~, mat_bin_idx] = histcounts(mat6090_use_pos, mat_edges);
[~, ~, map_bin_idx] = histcounts(map6090_use_pos, map_edges);

mat_change = mat9020_use_pos - mat6090_use_pos;
map_change = map9020_use_pos - map6090_use_pos;

% calculate the average arrow direction within each mat&map group
mean_mat_vec = zeros(num_bins, num_bins);
mean_map_vec = zeros(num_bins, num_bins);
counts = zeros(num_bins, num_bins);
for k = 1:numel(mat_bin_idx)
    i = mat_bin_idx(k);
    j = map_bin_idx(k);
    if i > 0 && j > 0 && i <= num_bins && j <= num_bins
        mean_mat_vec(i,j) = mean([mean_mat_vec(i,j), mat_change(k)], 'omitnan');
        mean_map_vec(i,j) = mean([mean_map_vec(i,j), map_change(k)], 'omitnan');
        counts(i,j) = counts(i,j) + 1;
    end
end
min_samples = 5;
mean_mat_vec(counts < min_samples) = NaN;
mean_map_vec(counts < min_samples) = NaN;

plot_points = (is_inside_squ & ~isnan(mean_map_vec) & ~isnan(mean_mat_vec));

% vector plot
q = quiver3(Px(plot_points), Ty(plot_points), 100*ones(size(Px(plot_points))), ...
    mean_map_vec(plot_points), ...
    mean_mat_vec(plot_points), ...
    zeros(size(mean_map_vec(plot_points))), ...
    1.6, ...
    'AutoScale','off', ...
    'Color','k',...
    'LineWidth',1.5);
    % 'MaxHeadSize',0.001,... 
q.ShowArrowHead = 'off';
hold on
xlim([0 4000]);
ylim([-5 35]);

% Calculate the conversion factor from data space to pixel space
ax = gca;
ax.Units = 'pixels';
ax_pos = ax.Position;
x_pixels_per_data = ax_pos(3)/diff(xlim);
y_pixels_per_data = ax_pos(4)/diff(ylim);

% triangle size
fixed_triangle_height_px = 2.8; % height (pixel)
fixed_triangle_width_px = fixed_triangle_height_px/5*4;  % width (pixel)

for i = 1:length(q.XData)
    x_start = q.XData(i);
    y_start = q.YData(i);
    x_end = q.XData(i) + q.UData(i)/1.15;
    y_end = q.YData(i) + q.VData(i)/1.15;
    
    % the vector
    v = [x_end - x_start, y_end - y_start];
    v_length = norm(v);

    if v_length == 0
        plot(x_end, y_end, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
        continue;
    end
    
    % unit vector
    v_unit = v / v_length;
    
    % perpendicular vector (based on pixel)
    v_px = [v(1)*x_pixels_per_data, v(2)*y_pixels_per_data];
    v_px_unit = v_px / norm(v_px);
    perp_px = [-v_px_unit(2), v_px_unit(1)]; 
    
    % vertices of a triangle (pixel)
    tip_px = [x_end*x_pixels_per_data, y_end*y_pixels_per_data];
    base_center_px = tip_px - v_px_unit * fixed_triangle_height_px * 0.8;
    base_left_px = base_center_px + perp_px * fixed_triangle_width_px/2;
    base_right_px = base_center_px - perp_px * fixed_triangle_width_px/2;
    
    % convert back to the data space
    tip = [tip_px(1)/x_pixels_per_data, tip_px(2)/y_pixels_per_data];
    base_left = [base_left_px(1)/x_pixels_per_data, base_left_px(2)/y_pixels_per_data];
    base_right = [base_right_px(1)/x_pixels_per_data, base_right_px(2)/y_pixels_per_data];
    
    % plot the triangles
    triangle_x = [tip(1), base_left(1), base_right(1), tip(1)];
    triangle_y = [tip(2), base_left(2), base_right(2), tip(2)];
    triangle_z = 100*ones(size(triangle_x));
    fill3(triangle_x, triangle_y, triangle_z, 'k', 'EdgeColor', 'k', 'LineWidth', 1);
end
set(gca, 'Fontname', 'Arial', 'Fontsize', 38);
xlabel('MAP (mm)', 'Fontsize', 41);
ylabel('MAT (^{o}C)', 'Fontsize', 41);
grid on;
box on;
set(gcf,'Position',[1 49 1920 954]);
ax = gca;
set(ax,'Position',[520,82,800,800]);


%% plot relative niche loss hotspot in MAT-MAP space 

reloss_tn =  mean(relativeloss_tmpniche,2,'omitnan') * (-100);  % relaitve niche loss
reloss_pn =  mean(relativeloss_preniche,2,'omitnan') * (-100);
reloss_tn_norm = normalize(reloss_tn,"range");
reloss_pn_norm = normalize(reloss_pn,"range");
nicheloss_tpcombined = (reloss_tn_norm + reloss_pn_norm) / 2;
[nicheloss_comb_lognorm,c,s] = normalize (log(nicheloss_tpcombined+1));

nicheloss_comb_lognorm = array2table([X,Y,nicheloss_comb_lognorm]);
nicheloss_comb_lognorm.Properties.VariableNames = {'Lon','Lat','nicheloss_comb_lognorm'};
writetable(nicheloss_comb_lognorm,'nicheloss_comb_lognorm_10km.csv')

normparam_10km = struct();
normparam_10km.mintn = min(reloss_tn);
normparam_10km.maxtn = max(reloss_tn);
normparam_10km.minpn = min(reloss_pn);
normparam_10km.maxpn = max(reloss_pn);
save normparam_10km normparam_10km

surf_smooth = 7;
heatmap_smooth = 7;
[Surf_niche_change] = niche2D_matmap(mat9020_use, map9020_use, nicheloss_tpcombined, 1, out_res, surf_smooth, heatmap_smooth); 
title('Relative niche loss','FontSize',10); 
ylim([-5 35]);

% a = array2table([mat9020_use, map9020_use, conflictfreq10km./pixelarea,nicheloss_tpcombined]);
% a.Properties.VariableNames = {'MAT','MAP','Conflict density','Relative niche loss'};
% maxRows = 900000; 
% part1 = a(1:maxRows, :);
% part2 = a(maxRows+1:2*maxRows, :);
% part3 = a(2*maxRows+1:end, :);
% writetable(part1, 'E:\DATA\Conflict\Manuscript\final version\source data\figure source data.xlsx', 'Sheet', 'Figure2_part1', 'WriteMode', 'overwritesheet');
% writetable(part2, 'E:\DATA\Conflict\Manuscript\final version\source data\figure source data.xlsx', 'Sheet', 'Figure2_part2', 'WriteMode', 'overwritesheet');
% writetable(part3, 'E:\DATA\Conflict\Manuscript\final version\source data\figure source data.xlsx', 'Sheet', 'Figure2_part3', 'WriteMode', 'overwritesheet');

hold on;
contour_color = colormap('Copper');
plot3(nichecontour_x(nichelevelgroup == 10),nichecontour_y(nichelevelgroup == 10),nichecontour_z(nichelevelgroup == 10),'Color',contour_color(1,:)  ,'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 9), nichecontour_y(nichelevelgroup == 9), nichecontour_z(nichelevelgroup == 9), 'Color',contour_color(10,:) ,'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 8), nichecontour_y(nichelevelgroup == 8), nichecontour_z(nichelevelgroup == 8), 'Color',contour_color(20,:) ,'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 7), nichecontour_y(nichelevelgroup == 7), nichecontour_z(nichelevelgroup == 7), 'Color',contour_color(60,:) ,'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 6), nichecontour_y(nichelevelgroup == 6), nichecontour_z(nichelevelgroup == 6), 'Color',contour_color(100,:),'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 5), nichecontour_y(nichelevelgroup == 5), nichecontour_z(nichelevelgroup == 5), 'Color',contour_color(140,:),'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 4), nichecontour_y(nichelevelgroup == 4), nichecontour_z(nichelevelgroup == 4), 'Color',contour_color(180,:),'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 3), nichecontour_y(nichelevelgroup == 3), nichecontour_z(nichelevelgroup == 3), 'Color',contour_color(200,:),'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 2), nichecontour_y(nichelevelgroup == 2), nichecontour_z(nichelevelgroup == 2), 'Color',contour_color(220,:),'LineWidth', 1.5);
plot3(nichecontour_x(nichelevelgroup == 1), nichecontour_y(nichelevelgroup == 1), nichecontour_z(nichelevelgroup == 1), 'Color',contour_color(250,:),'LineWidth', 1.5);

[~,niche_contour] = contour(outx,outy,Surf_niche_change');
maxnicheval = max(Surf_niche_change,[],'all');
contour_level_vec = [maxnicheval*0.9  maxnicheval*0.8  maxnicheval*0.7  maxnicheval*0.6];
niche_contour.LevelList = contour_level_vec;
mycolor = colormap('turbo');    
colormap(mycolor(15:end,:));
ylim([-5 35]); xlim([0 4000]);


%% Vector plot overlain on the niche loss heatmap

% mat map data
idx = (isnan(mat6090_use) | isnan(mat9020_use) | isnan(map6090_use) | isnan(map9020_use) | ...
    mat6090_use<mat_range(1) | mat6090_use>mat_range(2) | mat9020_use<mat_range(1) | mat9020_use>mat_range(2) |...
    map6090_use<map_range(1) | map6090_use>map_range(2) | map9020_use<map_range(1) | map9020_use>map_range(2));
loc_idx = logical(zeros(size(idx)));
loc_idx(conflictfreq10km_pos_loc) = true;
idx = (idx & ~loc_idx);
mat6090_use_pos = mat6090_use(~idx);
map6090_use_pos = map6090_use(~idx);
mat9020_use_pos = mat9020_use(~idx);
map9020_use_pos = map9020_use(~idx);

% calculate the arrow direction
[~, ~, mat_bin_idx] = histcounts(mat6090_use_pos, mat_edges);
[~, ~, map_bin_idx] = histcounts(map6090_use_pos, map_edges);

mat_change = mat9020_use_pos - mat6090_use_pos;
map_change = map9020_use_pos - map6090_use_pos;

% calculate the average arrow direction within each mat&map group
mean_mat_vec = zeros(num_bins, num_bins);
mean_map_vec = zeros(num_bins, num_bins);
counts = zeros(num_bins, num_bins);
for k = 1:numel(mat_bin_idx)
    i = mat_bin_idx(k);
    j = map_bin_idx(k);
    if i > 0 && j > 0 && i <= num_bins && j <= num_bins
        mean_mat_vec(i,j) = mean([mean_mat_vec(i,j), mat_change(k)], 'omitnan');
        mean_map_vec(i,j) = mean([mean_map_vec(i,j), map_change(k)], 'omitnan');
        counts(i,j) = counts(i,j) + 1;
    end
end
min_samples = 5;
mean_mat_vec(counts < min_samples) = NaN;
mean_map_vec(counts < min_samples) = NaN;

plot_points = (is_inside_squ & ~isnan(mean_map_vec) & ~isnan(mean_mat_vec));

% vector plot
q = quiver3(Px(plot_points), Ty(plot_points), 100*ones(size(Px(plot_points))), ...
    mean_map_vec(plot_points), ...
    mean_mat_vec(plot_points), ...
    zeros(size(mean_map_vec(plot_points))), ...
    1.6, ...
    'AutoScale','off', ...
    'Color','k',...
    'LineWidth',1.5);
    q.ShowArrowHead = 'off';
hold on
xlim([0 4000]);
ylim([-5 35]);

% Calculate the conversion factor from data space to pixel space
ax = gca;
ax.Units = 'pixels';
ax_pos = ax.Position;
x_pixels_per_data = ax_pos(3)/diff(xlim);
y_pixels_per_data = ax_pos(4)/diff(ylim);

% triangle size
fixed_triangle_height_px = 2.8; % height (pixel)
fixed_triangle_width_px = fixed_triangle_height_px/5*4;  % width (pixel)

for i = 1:length(q.XData)
    x_start = q.XData(i);
    y_start = q.YData(i);
    x_end = q.XData(i) + q.UData(i)/1.15;
    y_end = q.YData(i) + q.VData(i)/1.15;
    
    % the vector
    v = [x_end - x_start, y_end - y_start];
    v_length = norm(v);

    if v_length == 0
        plot(x_end, y_end, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
        continue;
    end
    
    % unit vector
    v_unit = v / v_length;
    
    % perpendicular vector (based on pixel)
    v_px = [v(1)*x_pixels_per_data, v(2)*y_pixels_per_data];
    v_px_unit = v_px / norm(v_px);
    perp_px = [-v_px_unit(2), v_px_unit(1)]; 
    
    % vertices of a triangle (pixel)
    tip_px = [x_end*x_pixels_per_data, y_end*y_pixels_per_data];
    base_center_px = tip_px - v_px_unit * fixed_triangle_height_px * 0.8;
    base_left_px = base_center_px + perp_px * fixed_triangle_width_px/2;
    base_right_px = base_center_px - perp_px * fixed_triangle_width_px/2;
    
    % convert back to the data space
    tip = [tip_px(1)/x_pixels_per_data, tip_px(2)/y_pixels_per_data];
    base_left = [base_left_px(1)/x_pixels_per_data, base_left_px(2)/y_pixels_per_data];
    base_right = [base_right_px(1)/x_pixels_per_data, base_right_px(2)/y_pixels_per_data];
    
    % plot the triangles
    triangle_x = [tip(1), base_left(1), base_right(1), tip(1)];
    triangle_y = [tip(2), base_left(2), base_right(2), tip(2)];
    triangle_z = 100*ones(size(triangle_x));
    fill3(triangle_x, triangle_y, triangle_z, 'k', 'EdgeColor', 'k', 'LineWidth', 1);
end
set(gca, 'Fontname', 'Arial', 'Fontsize', 38);
xlabel('MAP (mm)', 'Fontsize', 41);
ylabel('MAT (^{o}C)', 'Fontsize', 41);
grid on;
box on;
set(gcf,'Position',[1 49 1920 954]);
ax = gca;
set(ax,'Position',[520,82,800,800]);
title(''); 


