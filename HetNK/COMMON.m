function NewParams = COMMON(Params)
v2struct(Params);

[EMesh,BetaMesh,BMesh] = ndgrid(EGrid,BetaGrid,BGrid);
TauEMesh = ndgrid(TauEGrid,BetaGrid,BGrid);

EMesh = reshape(EMesh,[],BPts);
BetaMesh = reshape(BetaMesh,[],BPts);
BMesh = reshape(BMesh,[],BPts);
TauEMesh = reshape(TauEMesh,[],BPts);

[EDistMesh,BetaDistMesh,BDistMesh] = ndgrid(EGrid,BetaGrid,BDistGrid);
TauEDistMesh = ndgrid(TauEGrid,BetaGrid,BDistGrid);

EDistMesh = reshape(EDistMesh,[],BDistPts);
BetaDistMesh = reshape(BetaDistMesh,[],BDistPts);
BDistMesh = reshape(BDistMesh,[],BDistPts);
TauEDistMesh = reshape(TauEDistMesh,[],BDistPts);

sizeEVPrime = [EPts*BetaPts,BPts];
sizeDist = [EPts*BetaPts,BDistPts];

EBetaTrans = kron(BetaTrans,ETrans);

dimEV = repmat({':'},1,length(sizeEVPrime));
dimDist = repmat({':'},1,length(sizeDist));

clear Params;
NewParams = v2struct;
end