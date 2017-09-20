function aggRslt = aggregate_histogram(dist,w,r,params)
eGrid = params.eGrid;
kGrid = params.kGrid;

K = sum( reshape( reshape(kGrid,1,[]).*dist, [],1) );
L = sum( reshape( reshape(eGrid,[],1).*dist, [],1) );
aggRslt = v2struct(K,L);
end