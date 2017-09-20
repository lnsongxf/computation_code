function dY = solve_block(SL,S,SF,f)
n = size(S,1);
T = size(S,3);

dY = zeros(n,T);

% Compute C and d Matrix
C = zeros(n,n,T);
d = zeros(n,T);

% First period
C(:,:,1) = S(:,:,1)^(-1) * SF(:,:,1);
d(:,1) = -S(:,:,1)^(-1) * f(:,1);
for t=2:T
    C(:,:,t) = (S(:,:,t) - SL(:,:,t)*C(:,:,t-1))^(-1) * SF(:,:,t);
    d(:,t) = -(S(:,:,t)-SL(:,:,t)*C(:,:,t-1))^(-1) * (f(:,t)+SL(:,:,t)*d(:,t-1));
end

% Compute dY backwardly
dY(:,T) = d(:,T);
for t=T-1:-1:1
    dY(:,t) = d(:,t) - C(:,:,t)*dY(:,t+1);
end
end