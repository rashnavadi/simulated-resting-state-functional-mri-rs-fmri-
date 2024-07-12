function [static, dynamic] = hoco_fbglmfit(X, Obs, Wsize, param11, param12, param21, TR)
%   FBGLMFIT Fit a bayesian generalized linear model.
%   B = FBGLMFIT(X,OBS,WSIZE) fits a 2-level generalized linear model using the
%   predictor vector X to the response vector Obs. 
%   The parameters are estimated using Parametric Empirical Bayes.  

%   Noise distribution is assumed to be gaussian with zero mean and
%   non-spherity defined by Qs. 

%   Static and Dynamic Connectivity results are saved in outputs.  

%   X is a matrix with rows corresponding to observations, and columns to
%   predictor variables.

%   GLMFIT automatically includes a constant term in the
%   model by having shift equal to 1(do not enter a column of ones directly into X).  

%   Stats field contains 
%       'dfe'       degrees of freedom for error
%       's'         estimated dispersion parameter
%       'se'        standard errors of coefficient estimates B
%       'coeffcorr' correlation matrix for B
%       'covb'      estimated covariance matrix for B
%       't'         t statistics for B
%       'p'         p-values for B
%       'resid'     residuals

%   Example:  Fit a probit regression model for y on x.  Each y(i) is the
%   number of successes in n(i) trials.
%
%       x = [2100 2300 2500 2700 2900 3100 3300 3500 3700 3900 4100 4300]';
%       n = [48 42 31 34 31 21 23 23 21 16 17 21]';
%       y = [1 2 0 3 8 8 14 17 19 15 17 21]';
%       b = glmfit(x, [y n], 'binomial', 'link', 'probit');
%       yfit = glmval(b, x, 'probit', 'size', n);
%       plot(x, y./n, 'o', x, yfit./n, '-')
%
%   See also GLMFIT, GLMVAL, REGSTATS, REGRESS.

%   References:
%      [1] Friston, K.J. (2002) AClassical and Bayesian Inference in
%      Neuroimaging: Theory, Neuroimage 16, 465-483

%   Copyright 2014 Seaman Family MR Centre.
      
    if nargin < 3 
        Error(message('stats:gchi2cdf:TooFewInputs'));
    end 

    % you may want to check for correct X, Obs matrix size, Wsize later 
    npts = length(X);

    % Normalizing
    X = (X - mean(X)) ./ std(X);
    Obs = (Obs - mean(Obs)) ./ std(Obs); 
    

    % Static connectivity analysis  
    [static.bb,static.dev,static.stats] = glmfit(X,Obs,'normal');
    static.bb = atanh(static.bb); % Convert to Fisher's Z
    
    % Construct first-level design matrix 
    tmp_v = zeros(npts,1);
    tmp_v(1:(ceil(Wsize/2))) = 1; 
    DsgnMtx = toeplitz(tmp_v);
    DsgnMtx = DsgnMtx./(repmat(sum(DsgnMtx,2),1,npts)); % Vi from the paper
    X1 = DsgnMtx .* repmat(X,[1,npts]);


    %%%% Linear cov model for second level %%%%% 
    tmp        = zeros(1,npts);
    tmp(1:Wsize)  = 1/Wsize*(Wsize:-1:1);  % Linear covariance   
    tmp_covLINEAR = toeplitz(tmp); 

    % Set input params for spmn_PEB.m
    OPT = 1; 
    y = Obs; 
    P{1}.X = X1;
    P{2}.X = ones(npts,1);
    %P{2}.X = ones(4,1);
    P{1}.C{1} = param11 * eye(npts); % Model White noise  
    P{1}.C{2} = param12 * (toeplitz(exp(-(0:1:npts-1)*TR))); % to model AR(1) with exp(-TR)
    P{2}.C{1} = param21 * tmp_covLINEAR;

    clear tmp* 
    warning off; 

    %%%%%% Run spm.PEB  %%%%%%
    [C,P,F] = spm_PEB(y,P,OPT); 

    %%%%%% Prepare Output  %%%%%%
    
    dynamic.PEB.C = C;
    dynamic.PEB.P = P;
    dynamic.PEB.F = F;

    dynamic.bb = C{2}.E; 

    dynamic.stats.resid = Obs - (P{1}.X*C{2}.E); % residuals based on firsl level a posteriori information 
    dynamic.stats.beta = dynamic.bb;
    dynamic.stats.covb = C{2}.C;
    dynamic.stats.t = dynamic.stats.beta./sqrt(diag(dynamic.stats.covb));
    dynamic.stats.se = sqrt(diag(dynamic.stats.covb));
    dynamic.stats.t2 = C{3}.E/C{3}.C; % 2nd level beta value/SE second level value
    % Definitions down here need reconsideration 
    dynamic.dev = norm(dynamic.stats.resid)^2; 
    %dynamic.stats.dfe = static.stats.dfe; 
    %dynamic.stats.s = sqrt(dynamic.dev/dynamic.stats.dfe);    
    %dynamic.stats.covNse =  C{1}.M; 
    %dynamic.stats.se = sqrt(diag(dynamic.stats.covb));
    %dynamic.stats.coeffcorr = dynamic.stats.covb ./ (dynamic.stats.se*dynamic.stats.se');
     
    % Coefficient of Determination  
    dynamic.stats.R2 = (static.dev - dynamic.dev )/ static.dev;

    
    % Summary measures
    %Nsim = 1E5; 
    %rnd = mvnrnd(zeros(1,npts),single(full(dynamic.stats.covb)),Nsim); 
    %var_rnd = var(rnd');
    %var_con = var(dynamic.stats.beta); 
    
    %dynamic.summary.varPfa = numel(find(double(var_rnd)>var_con))/Nsim;
         
    dynamic.summary.connectivity = mean(full(dynamic.bb));
    dynamic.summary.conntStat    = dynamic.stats.t2;
    dynamic.summary.connSD      = std(full(dynamic.stats.beta));
    dynamic.summary.connRange    = range (full(dynamic.stats.t));    

    
%     % Variability Stats      
%     [tmp_rowRPT,tmp_colRPT] = meshgrid(mean(dynamic.stats.covb),mean(dynamic.stats.covb));
%     dynamic.summary.connvarCov = (dynamic.stats.covb - tmp_rowRPT - tmp_colRPT + mean(dynamic.stats.covb(:)))/(npts-1); % Unbiased estimation of variance, 1/N-1 factor
% 
%     [tmp_eigV,tmp_eigD] = eig(dynamic.summary.connvarCov);
%     tmp_eigd = diag(tmp_eigD);
%     tmp_Idx = find(cumsum(tmp_eigd)>0.80*sum(tmp_eigd));
%     tmp_eigd(tmp_Idx) = [];
%     dynamic.summary.NonSingCov = diag(tmp_eigd); % Non singular covariance matrix
%     dynamic.summary.npts = length(tmp_eigd);
%     dynamic.summary.varPvalue  = 1 - fgchi2cdf(zeros(dynamic.summary.npts,1) , dynamic.summary.NonSingCov, dynamic.summary.connvar, eye(dynamic.summary.npts), zeros(dynamic.summary.npts,1));
%     
%     clear tmp*       

end

