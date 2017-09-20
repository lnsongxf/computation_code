function Dist_new = UPDATE_DIST(Dist,bPolicyDist,Params)
BDistGrid = Params.BDistGrid;
BDistPts = Params.BDistPts;

EPts = Params.EPts;
BetaPts = Params.BetaPts;
EBetaTrans = Params.EBetaTrans;

sizeDist = Params.sizeDist;


% Look over grids
[~,bFutureLeftCell] = histc(bPolicyDist,BDistGrid);
% No extrapolation at the top
bFutureLeftCell = min(bFutureLeftCell,BDistPts-1);
bFutureRightCell = bFutureLeftCell + 1;
% Compute share
bFutureRightShare = (bPolicyDist-BDistGrid(bFutureLeftCell)) ./ (BDistGrid(bFutureRightCell)-BDistGrid(bFutureLeftCell));
bFutureLeftShare = 1-bFutureRightShare;

% Allocate Share to left and right
DistRight = zeros(sizeDist);
DistLeft = zeros(sizeDist);
for j=1:EPts*BetaPts
    DistRight(j,:) = accumarray(bFutureRightCell(j,:)', Dist(j,:)'.*bFutureRightShare(j,:)', [BDistPts 1]);
    DistLeft(j,:)  = accumarray(bFutureLeftCell(j,:)', Dist(j,:)'.*bFutureLeftShare(j,:)', [BDistPts 1]);
end

Dist_new = DistLeft + DistRight;
Dist_new = EBetaTrans'*Dist_new;
end