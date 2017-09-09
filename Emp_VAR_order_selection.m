%% Estimating Vector Autoregressions
% This task shows how to estimate a VAR(1) for three series.
% Source: https://www.kevinsheppard.com/MFE_MATLAB#Vector_Autoregressions
 
%% Prepare the workspace
% Clean up before starting
clear all
close all
clc
  
%% Estimate VAR(p) for = in 0,1,...,12
% The strategy is to estimate the VAR for all lag lengths between 0 and 12,
% and then to to compute the ICs and LRs.
%
% Note: the function vectorar was recently updated and so you might need a
% newer version to run this code.
%
% Note: The data length is adjusted to ensure that the same amount of data
% is used in all models.  This is important since ICs do not work correctly
% if the LHS data size is changing.

[data, xlstext] = xlsread('data_US');
s2 = cell(13,1);
AIC = zeros(13,1);
HQIC = zeros(13,1);
BIC = zeros(13,1);
LR = zeros(13,1);
pval = zeros(13,1);
for lags=0:12
    p = lags;
    disp('Data length after adjustment')
    T =length(data(12-lags+1:end,:)) - lags;
    disp(T)
    [~,~,~,~,~,~,~,~,s2{lags+1}] =  vectorar(data(12-lags+1:end,:),1,1:p);
    s2_det = det(s2{lags+1});
    AIC(lags+1) = log(det(s2_det)) + 3 + 2 * 3^2 * p /  T;
    HQIC(lags+1) = log(det(s2_det)) + 3 + 2 * log(log(T)) * 3^2 * p /  T;
    BIC(lags+1) = log(det(s2_det)) +  log(T) * 3^2 * p / T;
    if lags>=1
        LR(lags + 1) = (T - p * 3^2) * (log(det(s2{lags})) - log(det(s2{lags+1})));
        pval(lags+1) = 1 - chi2cdf(LR(lags+1),9);
    end
end
 
%% ICs
% For each IC, the minimum is selected.
disp('AIC selects')
disp(find(AIC==min(AIC)) - 1)
disp('HQIC selects')
disp(find(HQIC==min(HQIC)) - 1)
disp('BIC selects')
disp(find(BIC==min(BIC)) - 1)
 
%% Likelihood Ratio Testing
% The LR starts from a general model and selects down until the null is
% rejected.
disp([0:12;pval'])

