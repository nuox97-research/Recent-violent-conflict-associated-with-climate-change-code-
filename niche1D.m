function [outputArg1,outputArg2] = niche1D(inputArg1,inputArg2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% parameters for visualization
mat_range   = [-5.0 35.0];                              % data range of mean annual temperature 
matbin      = 2.0;                                      % moving window size along MAT (for running mean curves)
matmovbin   = 1.0;                                      % moving window step along MAT (for running mean curves)
matbinvec   = mat_range(1):matmovbin:mat_range(2);
matbinvec   = matbinvec + matmovbin/2.0;                % vector for drawing running mean curves
matvec      = -5:0.1:35; matvec = matvec.';             % vector for drawing fitted Gaussian curves
nichefactor = 1e6;                                      % scaling factor for normalizing niche (to avoid too low values)
leastsample = 600;

ft = fittype( 'gauss2' ); 
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off'; 
opts.Lower = [-Inf -Inf 0 -Inf 25 0]; 
opts.Upper = [Inf 25 Inf Inf 33 5];

% population data used                                  % add "& india<1" to exclude India monsoon region 
mat_use = inputArg1;
pop_use = inputArg2;

flt = pop_use>=0 & isnan(mat_use)==0;
mat_use = mat_use(flt);
pop_use = normsumunity(nichefactor,pop_use (flt));    % Population density data from the Gridded Population of the World (GPW v4) year 2015, normalize to sum unity
                 
% calculate running mean of population density along MAT gradient        
meanpdvec  = runningmeandata(pop_use, mat_use, matbin, matmovbin, mat_range,leastsample);
figure;
hold on;
plot(matbinvec, meanpdvec,'-r','LineWidth',2);
xlabel('MAT');ylabel('Running mean Pop Density');
% fit niche with double Gaussian
excl = isnan(meanpdvec)==0;   
[fitmatgaus, ~] = fit( matbinvec(excl).', meanpdvec(excl).', ft, opts );
niche_curve = fitmatgaus(matvec);
plot(matvec, niche_curve,'--b','LineWidth',2);
hold off;

outputArg1 = [matbinvec.' meanpdvec.'];
outputArg2 = fitmatgaus;

end

