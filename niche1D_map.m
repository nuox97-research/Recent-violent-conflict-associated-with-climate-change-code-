function [outputArg1,outputArg2] = niche1D_map(inputArg1,inputArg2)
%UNTITLED Summary of this function goes here

% parameters for visualization
                           
map_range   = [0 4000];                                 % data range of mean annual precipitation 
mapbin      = 200;                                      % moving window size along MAP (for running mean curves)
mapmovbin   = 100;                                      % moving window step along MAP (for running mean curves)
mapbinvec   = map_range(1):mapmovbin:map_range(2);
mapbinvec   = mapbinvec + mapmovbin/2.0;   
mapvec      = 0:1:4000; mapvec = mapvec.';              % vector for drawing fitted Gaussian curves
nichefactor = 1e6;                                      % scaling factor for normalizing niche (to avoid too low values)
leastsample = 100;

ft = fittype( 'gauss1' );  
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off'; 

% population data used                                  % add "& india<1" to exclude India monsoon region 
map_use = inputArg1;
pop_use = inputArg2;

flt = pop_use>=0 & isnan(map_use)==0;
map_use = map_use(flt);
pop_use = normsumunity(nichefactor,pop_use (flt));    % Population density data from the Gridded Population of the World (GPW v4) , normalize to sum unity
                 
% calculate running mean of population density along MAP gradient        
meanpdvec  = runningmeandata(pop_use, map_use, mapbin, mapmovbin, map_range,leastsample);
figure;
hold on;
plot(mapbinvec, meanpdvec,'-r','LineWidth',2);
xlabel('MAP');ylabel('Running mean Pop Density');

% fit niche with Gaussian
excl = isnan(meanpdvec)==0;   
[fitmapgaus, ~] = fit( log(mapbinvec(excl)+1).', meanpdvec(excl).', ft, opts );
niche_curve = fitmapgaus(log(mapvec+1));

plot(mapvec, niche_curve,'--b','LineWidth',2);
hold off;

outputArg1 = [mapbinvec.' meanpdvec.'];
outputArg2 = fitmapgaus;

end




