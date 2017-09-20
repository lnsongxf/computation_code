function [gini,sample] = gini_histogram(pdf,grid,n,seed)
% GINI_HISTOGRAM: sample from discrete pdf and compute gini
rng(seed);
un = rand(n,1);

[~,sampleIdx] = histc(un,[0 cumsum(pdf(:)')]);

sample = grid(sampleIdx);

gini = ginicoeff(sample);
end