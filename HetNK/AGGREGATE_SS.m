function [N,B,B2,tau,C,Y,Dividend,DistE,hours] = AGGREGATE_SS(bPolicyDist,nPolicyDist,cPolicyDist,Dist,Params,z,G,R,wage)
EDistMesh = Params.EDistMesh;
BDistMesh = Params.BDistMesh;
BDistPts = Params.BDistPts;
EPts = Params.EPts;

% Aggregate labor supply
N = sum(Dist(:).*nPolicyDist(:).*EDistMesh(:));
hours = sum(Dist(:).*nPolicyDist(:));

% Aggregate future asset
B = sum(Dist(:).*bPolicyDist(:));
B2 = sum(Dist(:).*BDistMesh(:));

% Taxes to balance government budget
DistE = sum(reshape(Dist,EPts,[]),2);
tau = (B + G - B/R) / DistE(3);

% Aggregate consumption
C = sum(Dist(:).*cPolicyDist(:));

% Aggregate output
Y = z*N;

% Dividend
Dividend = Y - wage*N;

%
end