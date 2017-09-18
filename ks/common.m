function newParams = common(params)
v2struct(params);
% Joint distribution
ZETrans = reshape(ZETrans, [ePts zPts ePts zPts]);
% permute ZETrans to (z,e) -> (z',e')
ZETrans = permute(ZETrans, [2 1 4 3]);
fullDiscountedTrans = zeros(zPts,ePts,betaPts,zPts,ePts,betaPts);
fullTrans = zeros(zPts,ePts,betaPts,zPts,ePts,betaPts);
for i_z=1:zPts
    for i_e = 1:ePts
        for i_beta = 1:betaPts
            for i_zp = 1:zPts
                for i_ep = 1:ePts
                    for i_betap = 1:betaPts
                        fullDiscountedTrans(i_z,i_e,i_beta,i_zp,i_ep,i_betap) = ...
                            betaGrid(i_betap)*ZETrans(i_z,i_e,i_zp,i_ep)*betaTrans(i_beta,i_betap);
                        fullTrans(i_z,i_e,i_beta,i_zp,i_ep,i_betap) = ZETrans(i_z,i_e,i_zp,i_ep)*betaTrans(i_beta,i_betap);
                    end
                end
            end
        end
    end
end
% FullDiscountedTrans = kron(BetaDiscountedTrans, ZETrans);
% % now transition matrix is in the order (z, e, beta) -> (z', e', beta')
fullDiscountedTrans = reshape(fullDiscountedTrans, [zPts*betaPts*ePts zPts*betaPts*ePts]);


% Conditional transition matrix (e,beta) -> (e',beta') | (z,z')
eBetaEPrimeBetaPrimeZZPrime = permute(fullTrans,[2,3,5,6,1,4]);
eBetaEPrimeBetaPrimeCondZZPrime = reshape(eBetaEPrimeBetaPrimeZZPrime,[ePts*betaPts,ePts*betaPts,zPts,zPts]);
eBetaEPrimeBetaPrimeCondZZPrime = eBetaEPrimeBetaPrimeCondZZPrime ./ sum(eBetaEPrimeBetaPrimeCondZZPrime,2);

clear params;
newParams = v2struct();
end
