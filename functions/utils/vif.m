function vifVals = vif(X, addIntercept)
% VIF  Compute variance–inflation factors for the columns of X.
%
%   vifVals = VIF(X)            returns a row vector of VIF values
%                               for each predictor in X.  X must
%                               already contain only the explanatory
%                               variables (no intercept column).
%
%   vifVals = VIF(X,true)       will first append an intercept column
%                               of ones to X, so that VIFs are computed
%                               for the original predictors after
%                               adjusting for the constant term.
%
%   --- Notes ----------------------------------------------------------
%   • The ith VIF is 1 / (1 – R_i²), where R_i² is the coefficient
%     of determination from regressing Xi on all other Xj ( j ≠ i ).
%   • A common rule of thumb flags VIF > 5 (or > 10) as problematic.
%   -------------------------------------------------------------------
%
% Yoondo (2025-05-06)

    if nargin < 2
        addIntercept = false;
    end
    if addIntercept
        X = [ones(size(X,1),1) X];
    end

    % Center columns to improve numeric stability (optional)
    X = X - mean(X,1);

    [nObs, nPred] = size(X);
    vifVals       = nan(1, nPred);

    for i = 1:nPred
        % Treat Xi as the response and the remaining predictors as regressors
        y_i   = X(:,i);
        X_oth = X(:,[1:i-1, i+1:end]);

        % Regress y_i on X_oth using QR factorization (numerically stable)
        beta     = X_oth \ y_i;              % same as:  beta = regress(y_i, X_oth);
        y_hat    = X_oth * beta;
        SSE      = sum((y_i - y_hat).^2);
        SST      = sum((y_i - mean(y_i)).^2);
        R2_i     = 1 - SSE / SST;

        vifVals(i) = 1 / (1 - R2_i);
    end
end
