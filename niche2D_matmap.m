function [outputArg1] = niche2D_matmap(inputArg1,inputArg2, inputArg3, inputArg4,inputArg5,inputArg6,inputArg7)

matuse  = inputArg1;                                  % current mean annual precipitation (MAP) data 
mapuse  = inputArg2;                                  % current mean annual temperature (MAT) data
zfactor = inputArg3;                                  % density (population, conflict etc.)
outputsurface = inputArg4;                            % Y/N output surface
out_res = inputArg5;                                  % resolution of canvas visualizing heat maps 

smth = inputArg6;           % smooth output surface window size:          5/7 as default
heat_neighbor = inputArg7;  % heat map neighbor window size: 5/7 as default

flt = ~isnan(matuse) & ~isnan(mapuse) & ~isnan(zfactor);

map_range = [0 4000];                               % data range mean annual temperature
mat_range = [-10 40];                               % data range mean annual precipitation
datarange =  [map_range(1) map_range(2) mat_range(1) mat_range(2)];

% leastsample = 1000;                                 % a threshold for excluding sparse data points
CI_defined = 0.99;                                  % bins with point densities lower than 1% of maximum are excluded 

xx = linspace(map_range(1),map_range(2),out_res);                 
yy = linspace(mat_range(1),mat_range(2),out_res);
[outx, outy] = meshgrid(xx,yy);

% point density on a log-scale, producing Fig 1G
figure;                 
[~,pdensity] = cloudPlot(mapuse(flt),matuse(flt),datarange,true,[out_res,out_res]); 
ylabel ('MAT (^{o}C)','Fontsize',16); xlabel ('MAP (mm)','Fontsize',16); title ('Avaialbe niche (log point density)'); 
colormap('jet'); axis square; set(gca,'FontSize',16);

% remove bins with sparse data points to account for potential bias
pdensity = (pdensity'); pdensity( isnan(pdensity) ) = 0;
smoothsize = [smth,smth];     % smooth 5 by 5
Hfilter = fspecial('average',smoothsize);
pdensity = imfilter(pdensity,Hfilter); 
ci_flt = pdensity<max(pdensity(:))*(1.0-CI_defined); 
xdata = [map_range(1) map_range(2)]; ydata = [ mat_range(1) mat_range(2)];
Pointvalue = impixel(xdata,ydata,ci_flt,mapuse(flt),matuse(flt)); 
removepoint = Pointvalue == 1;
removepoint = removepoint(:,1);
px = mapuse(flt); py = matuse(flt); pz = zfactor(flt);
px(removepoint) = []; py(removepoint) = []; pz(removepoint) = []; % excluding data with too low point densities

% heat map of focal variable  with low density data excluded
[~,canvas_heatmap] = heatmapcompute(px,py,pz,datarange,false,false,[out_res,out_res],true,[heat_neighbor,heat_neighbor]); % compute mean population density each bin within the MAP-MAT space
colormap('jet');
canvas_heatmap = canvas_heatmap'; 
canvas_heatmap( isnan(canvas_heatmap) ) = 0;
canvas_heatmap = imfilter(canvas_heatmap,Hfilter);

% final smoothed heat map of focal variable
if outputsurface
    figure;
    surf(outx,outy,canvas_heatmap);
    colormap('jet'); axis square; set(gca,'FontSize',16);shading interp; view(2);
    ylabel ('MAT (^{o}C)','Fontsize',16); xlabel ('MAP (mm)','Fontsize',16); title ('Heat map density');
    xlim(map_range); ylim(mat_range);
end

outputArg1 = canvas_heatmap';

end

