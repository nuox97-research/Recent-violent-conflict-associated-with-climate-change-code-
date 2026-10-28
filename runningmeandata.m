%% a function used for calculating the running mean of population density along the gradient of mean annual temperature

function y = runningmeandata(pop_use, mat_use, bin, movbin, mat_range,leastsample)

start = mat_range(1);       
final = mat_range(2);        
binvec = start:movbin:final;
binvec = binvec+movbin/2;
iterator=start:movbin:final;
length_i=length(iterator);

meanpdvec  = zeros(1,length_i);
for i = 1:length_i
    j=iterator(i);
    databin = pop_use (mat_use>=j & mat_use<j+bin);
    n = length(databin);     
    if n<leastsample
        mean_pd = nan;
    else
        mean_pd = mean (databin,'omitnan');
    end
    meanpdvec(i) = mean_pd;
end

y = meanpdvec;


