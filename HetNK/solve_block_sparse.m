function J = solve_block_sparse(SL,S,SF,f)
n = size(S,1);
T = size(S,3);

a = cell(1,T);
for t=1:T
    a{t} = sparse(S(:,:,t));
end
A = blkdiag(a{:});

% +1 Off diagnoal
b = cell(1,T-1);
for t=1:T-1
    b{t} = sparse(SF(:,:,t));
end
B = blkdiag(b{:});

% -1 Off diagonal
c = cell(1,T-1);
for t=1:T-1
    c{t} = sparse(SL(:,:,t+1));
end
C = blkdiag(c{:});


% Stack the large matrix
J = A + [sparse(n*(T-1),n),B;sparse(n,n*T)] + [sparse(n,n*T);C,sparse(n*(T-1),n)];
end