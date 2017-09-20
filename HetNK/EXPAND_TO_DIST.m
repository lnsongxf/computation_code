function [bPolicyDist,nPolicyDist,cPolicyDist] = EXPAND_TO_DIST(bPolicy,nPolicy,wage,R,tau,Dividend,Params)
BGrid = Params.BGrid;
BPts = Params.BPts;
BDistGrid = Params.BDistGrid;
sizeDist = Params.sizeDist;
BDistMesh = Params.BDistMesh;
TauEDistMesh = Params.TauEDistMesh;
EDistMesh = Params.EDistMesh;

% Expand to policy to finer grids
bPolicyPp = pchip(BGrid,bPolicy);
bPolicyDist = myppual(bPolicyPp,BDistGrid);

nPolicyPp = pchip(BGrid,nPolicy);
nPolicyDist = myppual(nPolicyPp,BDistGrid);

cPolicyDist = BDistMesh - tau*TauEDistMesh + Dividend + wage*EDistMesh.*nPolicyDist - bPolicyDist/R;

end