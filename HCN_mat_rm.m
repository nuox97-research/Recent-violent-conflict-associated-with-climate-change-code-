function  HCN_mat_rm(inputArg1,inputArg2,inputArg3,inputArg4,inputArg5,inputArg6,inputArg7)

zfactor     = inputArg1;   
matuse      = inputArg2; 
mat_range   = inputArg3; 
leastsample = inputArg4;
ncolor      = inputArg5;
scalefactor = inputArg6;
ifnormsumunity = inputArg7;

nichefactor = 1e6;
[~,nmatuse] = size(matuse);
HCN_timeperiod = [];
for imat = 1:1:nmatuse
    iimatuse = matuse(:,imat);
    flt = ~isnan(iimatuse);
   
    if ifnormsumunity
        zfac_use = normsumunity(nichefactor,zfactor (flt));
    else
        zfac_use = zfactor (flt);
    end

    mean_zfac_vec  = runningmeandata(zfac_use, iimatuse(flt), 2.0, 1.0, mat_range,leastsample);
    HCN_timeperiod = [HCN_timeperiod; mean_zfac_vec];
end
maxcur_out = prctile(HCN_timeperiod,95)           *scalefactor;                  % 95 percentile of running mean ensemble (here equivalent to maximum due to sparse data points)
mincur_out = prctile(HCN_timeperiod,5)            *scalefactor;                   %  5 percentile of running mean ensemble (here equivalent to minimum due to sparse data points)
meancurRM_out = mean(HCN_timeperiod,'omitnan')    *scalefactor;                     %    mean value of running mean ensemble

%%% plotting Figure 1a
matbin      = 2.0;                                      % moving window size along MAT (for running mean curves)
matmovbin   = 1.0;                                      % moving window step along MAT (for running mean curves)
matbinvec   = mat_range(1):matmovbin:mat_range(2);
matbinvec   = matbinvec + matmovbin/2.0;   
% ft = fittype( 'gauss2' );  
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off'; opts.Lower = [-Inf -Inf 0 -Inf 25 0]; opts.Upper = [Inf Inf Inf Inf 30 5];        
matvec      = -5:0.1:35; 
matvec = matvec.';  
binvec = mat_range(1):1.0:mat_range(2)+0.5;
% figure;
% h_RM = gcf;
% hold on;
C = linspecer(9,'sequential'); 


fltout = isnan(meancurRM_out)==0;
h = fill([binvec(fltout),fliplr(binvec(fltout))],[mincur_out(fltout),fliplr(maxcur_out(fltout))],C(ncolor,:));
set(h,'LineStyle','none','facealpha',0.3);
hold on;
plot(binvec,meancurRM_out,'-','LineWidth',2, 'Color',C(ncolor,:));
xlabel ('MAT (^{o}C)','Fontsize',20);
% ylabel ('Population density','Fontsize',20); 
xlim([-5 35]); 
% ylim([0 1.75]);
box on;
set(gca,'YTickLabel',[],'ytick',[],'FontSize',16);
% title('Human niche','Fontsize',12);

end