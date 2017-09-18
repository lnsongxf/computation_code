function newParams = common(params)
v2struct(params);
% discount transition matrix of beta
% betaDiscountedTrans = betaTrans .* repmat(betaGrid(:), 1, betaPts);
% betaDiscountedTrans = betaTrans .* repmat(betaGrid(:)', betaPts, 1);

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

% conditional transition matrix e -> e' | z
eEPrimeZZPrime = permute(reshape(ZETrans, [zPts ePts zPts ePts]), [2 4 1 3]);
eTransConditionalOnZ = reshape(sum(eEPrimeZZPrime, 4), [ePts ePts zPts]);

% Conditional transition matrix (e,beta) -> (e',beta') | z
eBetaEPrimeBetaPrimeZZprime = permute(fullTrans,[2,3,5,6,1,4]);
eBetaTransConditionalOnZ = reshape(sum(eBetaEPrimeBetaPrimeZZprime,6),[ePts*betaPts,ePts*betaPts,zPts]);

% unconditional transition matrix
eBetaTrans = sum(eBetaTransConditionalOnZ,3)/2;
eBetaTransInv = eBetaTrans^1000;
eBetaTransInv = eBetaTransInv(1,:);

clear Params;
newParams = v2struct();
end
