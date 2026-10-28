function y = normsumunity(nichefactor, data_pop)

pop_sum = nansum(data_pop);  
pop_sumunity = nichefactor * data_pop./pop_sum; %% normalized to sum unity

y = pop_sumunity;