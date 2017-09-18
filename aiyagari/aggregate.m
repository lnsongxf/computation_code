function aggRslt = aggregate(dist,rn,w,r,params)
eTrans = params.eTrans;
eGrid = params.eGrid;

K_t = mean(dist.k_t);
eTransStationary = eTrans^1000;
eTransStationary = eTransStationary(1,:);

L = sum(eTransStationary.*eGrid);

aggRslt = v2struct(K_t,L);
end