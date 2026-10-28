function HCN_map(inputArg1,inputArg2,inputArg3,inputArg4,inputArg5,inputArg6,inputArg7)

zfactor = inputArg1;   
mapuse  = inputArg2; 
map_range   = inputArg3; 
leastsample = inputArg4;
ncolor = inputArg5;
scalefactor = inputArg6;
ifnormsumunity = inputArg7;

nichefactor = 1e6;

[~,nmapuse] = size(mapuse);
HCN_timeperiod = [];
for imap = 1:1:nmapuse
    iimapuse = mapuse(:,imap);
    flt = ~isnan(iimapuse);
    if ifnormsumunity
        zfac_use = normsumunity(nichefactor,zfactor (flt));
    else
        zfac_use = zfactor (flt);
    end

    mean_zfac_vec  = runningmeandata(zfac_use, iimapuse(flt), 500, 100, map_range,leastsample);
    HCN_timeperiod = [HCN_timeperiod; mean_zfac_vec];
end
maxcur_out = prctile(HCN_timeperiod,95) * scalefactor;                  % 95 percentile of running mean ensemble (here equivalent to maximum due to sparse data points)
mincur_out = prctile(HCN_timeperiod,5)  * scalefactor;                   %  5 percentile of running mean ensemble (here equivalent to minimum due to sparse data points)
meancurRM = mean(HCN_timeperiod,'omitnan')   * scalefactor;                     %    mean value of running mean ensemble

%%% plotting 
mapbin      = 500;                                      % moving window size along MAT (for running mean curves)
mapmovbin   = 100;                                      % moving window step along MAT (for running mean curves)
mapbinvec   = map_range(1):mapmovbin:map_range(2);
mapbinvec   = mapbinvec + mapmovbin/2.0;   
       
% mapvec      = 0:100:4000; 
% mapvec = mapvec.';  
binvec = map_range(1):100:map_range(2) + 50;
C = linspecer(9,'sequential'); 

% 
fltout = isnan(meancurRM)==0;
h = fill([binvec(fltout),fliplr(binvec(fltout))],[mincur_out(fltout),fliplr(maxcur_out(fltout))],C(ncolor,:));
set(h,'LineStyle','none','facealpha',0.3);
hold on;
plot(binvec,meancurRM,'-','LineWidth',2, 'Color',C(ncolor,:));

xlabel ('MAP (mm)','Fontsize',20);
% ylabel ('Population density','Fontsize',20); 
xlim([0 4000]); 
% ylim([0 1.75]);
box on;
set(gca,'YTickLabel',[],'ytick',[],'FontSize',16);
% title('Human niche precipitaiton','Fontsize',12);

end