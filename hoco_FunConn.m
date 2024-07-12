function varargout = hoco_FunConn(Time_series_file, WSize, SurrOption, ...
    SaveSuffix)
% Calculate the adjacency matrix using the time-series file provided by the
% user and Pearson's correlation coefficient.

% University of Calgary, 2019
% Raphael Fernandes Casseb
% =========================================================================
% INPUTS
% ------
% Time_series_file (string) Path to the file containing the time--series.
%                           File must be a 2xn cell array. First row 
%                           contains time-series IDs and second row stores 
%                           the time-series. One time series per column. 
%                           Time grows downwards.
%                           Eg:
%                               |   1   |   2   |   3   |   4   |   5   |
%                           ---------------------------------------------
%                           1   |'PCC'  |'MPFC' |'R_hip'|'L_hip'|'L_par'|
%                           ---------------------------------------------
%                           2   |180x1  |180x1  |180x1  |180x1  |180x1  |
%
%                           180x1 is a column vector of 180 time-points.

% SurrOption       (double) 0: Use no surrogate data methods. It is the 
%                              same as doing the standard cross correlation
%                              analysis.
%                           1: Create surrogate data and run the analysis 
%                              on them. 
%                              Surrogate data is created by Resampling both 
%                              time-series to create another two, such that
%                              no pair of data points is the same as the 
%                              original. This is done for the whole time-
%                              series (static FC) and also window-wise 
%                              (dynamic FC).

% SaveSuffix       (string) Suffix added to the saved files 


% OUTPUTS
% -------
% varargout{1}     (double) EffectiveWSz: the effective window size. This
%                           is important for when tapering is used (which
%                           is the default in this algorithm). The
%                           effective window size is always even when
%                           tapering is applied. 
%                           Eg: a user defined window size of 45 will yield
%                           an effective window size of 46, using WAlpha=3.

% =========================================================================
fprintf('Starting Pearson''s FC analysis: %s\n', datetime('now'))
[~, FileName] = fileparts(fileparts(Time_series_file));
fprintf('File name: %s\n', FileName)

% Additional inputs
% -----------------
is_symm = 1;  % 1 means the adjacency matrix is symmetric, so we can do
              % calculations only for the upper triangle (since the lower
              % is symmetric).
Taper  = 1;   % Taper window template? 1: yes; 0: no. Default is 1.
WAlpha = 3;   % Window alpha: inclination factor of the window edges. The 
              % smaller the alpha, the more rectangular the window.


% -------------------------------------------------------------------------
% 1. Prepare table to receive the adjancency matrix
% -------------------------------------------------------------------------
% Load time-series file
load(Time_series_file);
TSsData = TSs; clear TSs
nTSs = size(TSsData,2); % number of time-series
nTPs = size(TSsData{2,1},1); % number of time-points

% Initialize table with some info
% -------------------------------
StaticFCArray  = cell(nTSs+1,nTSs+1);
SldWFCArray    = cell(nTSs+1,nTSs+1);
SldWFC_SDArray = cell(nTSs+1,nTSs+1);

% Row and Column names
StaticFCArray(2:end,1)  = TSsData(1,:)';
SldWFCArray(2:end,1)    = TSsData(1,:)';
SldWFC_SDArray(2:end,1) = TSsData(1,:)';

StaticFCArray(1,2:end)  = TSsData(1,:); 
SldWFCArray(1,2:end)    = TSsData(1,:);
SldWFC_SDArray(1,2:end) = TSsData(1,:);

% Save subject identifier in the matrices
[Path, Sbj] = fileparts(fileparts(Time_series_file));
StaticFCArray{1,1}  = Sbj;
SldWFCArray{1,1}    = Sbj;
SldWFC_SDArray{1,1} = Sbj;


% -------------------------------------------------------------------------
% 2. Calculate connectivity matrices
% -------------------------------------------------------------------------
% 2.1. Initialize window template
if Taper
    WTemplate = GetWindowTemplate(nTPs, WAlpha, WSize);  % Creates gaussian sliding window with window alpha
    % Zero values below half maximum to avoid long roll offs
    WTemplate(WTemplate<0.5) = NaN;
    nWindows = nTPs - WSize + 1;
    if mod(WSize,2) % Odd 
        nWindows = nWindows - 1;
    end
    WindowMatrixTemplate = NaN(nTPs, nWindows);
    for i=1:nWindows
        WShift = circshift(WTemplate, round(-nTPs/2) + round(WSize/2) + ...
            (i-1));
        WindowMatrixTemplate(:,i) = WShift;
    end
    
else
    WTemplate = [ones(WSize,1); NaN(nTPs-WSize,1)];
    nWindows = nTPs - WSize + 1;
    WindowMatrixTemplate = NaN(nTPs, nWindows);
    for i=1:nWindows
        WindowMatrixTemplate(:,i) = circshift(WTemplate, i-1);
    end
    
end


EffectiveWSz = sum(~isnan(WTemplate));
EffectiveNWindows = size(WindowMatrixTemplate,2);
varargout{1} = EffectiveWSz;

% Loop through all time-series
for ii=1:nTSs 
    
    fprintf('....Working on ROI %d/%d\n', ii, nTSs)
    
    % Get 1st ROI TS
    X = TSsData{2,ii};
    X = (X-mean(X))/std(X);
    if SurrOption == 1
        [XMat, Rspd_WndwTmplMtx2] = ...
            ResampleMtxOfWindowTemplates(EffectiveWSz, ...
            EffectiveNWindows, WindowMatrixTemplate, X); 
    else
        XMat = WindowMatrixTemplate .* repmat(X,[1,nWindows]);
    end
    

    % Run one TS against all the other TSs
    % ------------------------------------
    for iii = 1:nTSs
        
        if ii==iii % Exclude diagonal elements
            continue
            
        else
            if is_symm % Exclude upper triangle
                if ii < iii
                    continue
                end
            end
            
            % fprintf('%d, %d\n',ii,iii)
            % Get 2nd ROI ts (modeled)
            Obs    = TSsData{2,iii};
            Obs    = (Obs-mean(Obs))/std(Obs);
            ObsMat = zeros(nTPs,EffectiveNWindows);
            if SurrOption == 1
                for i = 1:EffectiveNWindows
                    ObsMat(i:(i+EffectiveWSz-1),i) = Obs(Rspd_WndwTmplMtx2(:,i));
                end
                ObsMat = WindowMatrixTemplate.*ObsMat;
                
            else
                ObsMat = WindowMatrixTemplate.*...
                    repmat(Obs,1,EffectiveNWindows);
            end
   
            % Pearson's dynamic correlation
            %           -------
            SldWFCArray{ii+1,iii+1}    = arrayfun(@(k) atanh(corr(...
                XMat( ~(isnan(XMat(:,k))) ,k), ...
                ObsMat( ~(isnan(ObsMat(:,k))) ,k))), ...
                1:nWindows,'Uni',1)';
            SldWFC_SDArray{ii+1,iii+1} = std(SldWFCArray{ii+1,iii+1});
            
            % Pearson's static correlation
            %           ------
            if SurrOption == 1
                [X,~] = ResampleTimepoints(X,[]);
            end
            StaticFCArray{ii+1,iii+1} = atanh(corr(X,Obs));
            
        end
        
    end
    
end


% -------------------------------------------------------------------------
% 3. Save matrices
% -------------------------------------------------------------------------
save(fullfile(Path,Sbj,['mtx_StaticFC' SaveSuffix '.mat']),'StaticFCArray')
save(fullfile(Path,Sbj,['mtx_SldWFC' SaveSuffix '.mat']),'SldWFCArray')
save(fullfile(Path,Sbj,['mtx_SldWFC_SD' SaveSuffix '.mat']),'SldWFC_SDArray')
% toc


fprintf('Finished Pearson''s FC analysis: %s\n', datetime('now'))

end
% =========================================================================
% SUPPORT FUNCTIONS
% =========================================================================
% -------------------------------------------------------------------------
% A. Create Tapering Template
% -------------------------------------------------------------------------
function TaperTemplate = GetWindowTemplate(nTPs, WAlpha, WSize)

nT1 = nTPs;
if mod(nTPs, 2) ~= 0
    nTPs = nTPs + 1;
end

m = nTPs/2;
w = round(WSize/2);

x=0:nTPs-1;
gw = exp(- ((x-m).^2)/ (2 * WAlpha * WAlpha))';
b = zeros(nTPs, 1);  b((m -w + 1):(m+w)) = 1;
TaperTemplate = conv(gw, b); 
TaperTemplate = TaperTemplate/max(TaperTemplate); 
TaperTemplate = TaperTemplate(m+1:end-m+1);
TaperTemplate = TaperTemplate(1:nT1);

end

% -------------------------------------------------------------------------
% B. Resample All Time-Series Points with Replacement
% -------------------------------------------------------------------------
function [RspdTS1,RspdTS2] = ResampleTimepoints(TS1,TS2)

if isempty(TS2)
    Idx     = [1:size(TS1,1)]';
    IdxTS1   = [Idx TS1];
    ShfdIdx = datasample(IdxTS1,size(IdxTS1,1));
    % Check if any of the ShfdIdx vector remained at its original position
    while any(ShfdIdx(:,1) == Idx)
        IdxRep = find(ShfdIdx(:,1) == Idx);
        for ii = 1:size(IdxRep,1)
            ShfdIdx(IdxRep(ii),:) = datasample(IdxTS1,1);
        end
    end
    RspdTS1 = ShfdIdx(:,2);
    RspdTS2 = [];
    
else
    Idx     = [1:size(TS1,1)]';
    IdxTS1   = [Idx TS1];
    IdxTS2   = [Idx TS2];
    ShfdIdx1 = datasample(IdxTS1,size(IdxTS1,1));
    ShfdIdx2 = datasample(IdxTS2,size(IdxTS2,1));
    % Check if any of the ShfdIdx vector remained at its original position
    while any(ShfdIdx1(:,1) == ShfdIdx2(:,1))
        IdxRep = find(ShfdIdx1(:,1) == ShfdIdx2(:,1));
        for ii = 1:size(IdxRep,1)
            ShfdIdx1(IdxRep(ii),:) = datasample(IdxTS1,1);
        end
    end
    RspdTS1 = ShfdIdx1(:,2);
    RspdTS2 = ShfdIdx2(:,2);
end
    
end

% -------------------------------------------------------------------------
% C. Resample Time-Series Points with Replacement Within the Matrix of
% Window Templates
% -------------------------------------------------------------------------
function [Rspd_WndwMtx4TS, Rspd_WndwTmplMtx2] = ...
    ResampleMtxOfWindowTemplates(EffectiveWSz, EffectiveNWindows, ...
    WindowMatrixTemplate, TS)

        cc_Idx = ndgrid(1:EffectiveNWindows,1:EffectiveWSz)';
    
        % Full template of the window indexes
        WndwTmplMtx_Idx = repmat([1:EffectiveWSz]',[1 EffectiveNWindows]) ...
            + cc_Idx-1;
        
        % Resample window indexes with repetition
        Rspd_WndwTmplMtx1 = WndwTmplMtx_Idx(...
            randi(size(WndwTmplMtx_Idx,1),EffectiveWSz,...
            size(WndwTmplMtx_Idx,2)) + ...
            (0:size(WndwTmplMtx_Idx,2)-1)*size(WndwTmplMtx_Idx,1));
        Rspd_WndwTmplMtx2 = WndwTmplMtx_Idx(...
            randi(size(WndwTmplMtx_Idx,1),EffectiveWSz,...
            size(WndwTmplMtx_Idx,2)) + ...
            (0:size(WndwTmplMtx_Idx,2)-1)*size(WndwTmplMtx_Idx,1));
        
        % Keep resampling if there are equally positioned elements
        while any(any(Rspd_WndwTmplMtx1 == Rspd_WndwTmplMtx2))
            IdxRep = find(Rspd_WndwTmplMtx1 == Rspd_WndwTmplMtx2);
            [IdxRepR,IdxRepC] = ind2sub(size(Rspd_WndwTmplMtx1),IdxRep);
            for i=1:size(IdxRepR)
                Rspd_WndwTmplMtx1(IdxRepR(i),IdxRepC(i)) = ...
                    WndwTmplMtx_Idx(randi(EffectiveWSz),IdxRepC(i));
            end
        end
        
        % Create XMat
        Rspd_WndwMtx4TS = zeros(size(TS,1),EffectiveNWindows);
        for i = 1:EffectiveNWindows
            Rspd_WndwMtx4TS(i:(i+EffectiveWSz-1),i) = TS(Rspd_WndwTmplMtx1(:,i));
        end
        Rspd_WndwMtx4TS = Rspd_WndwMtx4TS.*WindowMatrixTemplate;
        
        % The same as the for above, but slower
        % TSMat = repmat(TS,1,EffectiveNWindows);
        % ind   = sub2ind(size(TSMat),FullWTemplate_RspIdx1,cc_Idx);
        % Rspd_WndowMtx4TS = TSMat(ind);
        
end