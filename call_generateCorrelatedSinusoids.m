clear;
close all;
% Initialize variables
% Number of iterations
numIterations = 10;

Basefile= '/Users/trashnavadi/Documents/Data_Analysis/2022/analyses/kmeans_investigation/2024/June/FINAL';

nSubj = numIterations;
% Timepoints
nTimePts = 150;
nStates = 2;
k = nStates;
TR = 2; % in sec

% WSize = [30, 40, 50]; % in timepoints (TR = 2sec)
WSize = 50; % in timepoints (TRs), 300 sec
window_size = WSize;

% transitionPoint = [15, 45, 75, 120] % in TR, which are [30, 90, 150, 240]
transitionPoint = 120;  % Transition between states happens at this point

% Preallocate arrays to store correlations
StateFlip_SWC = cell(nSubj,1);
StateFlip_HOCo = cell(nSubj,1);

all_corr_xy = zeros(1, numIterations);
all_corr_xz = zeros(1, numIterations);
all_corr_yz = zeros(1, numIterations);

all_corr_over_time_xy = zeros((nTimePts - WSize + 1), numIterations);
all_corr_over_time_xz = zeros((nTimePts - WSize + 1), numIterations);
all_corr_over_time_yz = zeros((nTimePts - WSize + 1), numIterations);

all_x = cell(1, numIterations);
all_y = cell(1, numIterations);
all_z = cell(1, numIterations);


% Run the function for 1000 iterations
for iter = 1:numIterations
%    [x, y, z, corr_over_time_xy, corr_over_time_xz, corr_over_time_yz] = generateCorrelatedSinusoids(nTimePts, transitionPoint);
   [x, y, z, corr_xy, corr_xz, corr_yz, corr_over_time_xy, corr_over_time_xz, corr_over_time_yz] = generateCorrelatedSinusoids(nTimePts, transitionPoint);

   % Save x, y, and z for each iteration
   all_x{iter} = x;
   all_y{iter} = y;
   all_z{iter} = z;
   % Store the correlation values for this iteration
   all_corr_xy(iter) = corr_xy;
   all_corr_xz(iter) = corr_xz;
   all_corr_yz(iter) = corr_yz;

   all_corr_over_time_xy(:,iter) = corr_over_time_xy';
   all_corr_over_time_xz(:,iter) = corr_over_time_xz';
   all_corr_over_time_yz(:,iter) = corr_over_time_yz';


   clear TSs;TSs = cell(2,3);
   TSs{1,1} = 'timeseires_1';
   TSs{2,1} = x;
   TSs{1,2} = 'timeseires_2';
   TSs{2,2} = y;
   TSs{1,3} = 'timeseires_3';
   TSs{2,3} = z;
   save(fullfile(Basefile, '/sim_timeseries_1.mat'), 'TSs')

   load(fullfile(Basefile,'/sim_timeseries_1.mat'))
   SaveSuffix='Pearsons'; SurrOption=0;
   Time_series_file = fullfile(Basefile,'/sim_timeseries_1.mat');
   hoco_FunConn(Time_series_file, WSize, SurrOption, SaveSuffix)
   fprintf('calculating the Pearsons FC matrices between each pairs of the three concatenated timeseries \n')

   % the output from hoco_FunConn will be saved in the dir
   load(fullfile(Basefile, 'mtx_SldWFCPearsons.mat'))
   FC_21 = SldWFCArray{3,2};
   FC_31 = SldWFCArray{4,2};
   FC_32 = SldWFCArray{4,3};
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Step 6: Apply k-means clustering
   k = 2;  % Number of clusters
   % kmeans on the full concatenated timeseries
   kmeans_data_SWC = [FC_21, FC_31, FC_32];
   [StateIdx_SWC, Centroid_SWC, SumDistPartial_SWC] = kmeans(kmeans_data_SWC, k, 'Distance', 'cityblock', 'MaxIter',100, 'Replicates', 100);

   % record the time of transition between states
   % -------------
   StateFlip_SWC{iter} = find(logical(diff(StateIdx_SWC)));
   % ========================================
   X_full   = TSs{2,1};
   Obs_full = TSs{2,2};
   T_full   = TSs{2,3};
   % run HOCo on the concatenated timeseries
   param11 = 1; param12 = 1; param21 = 1;

   % HOCo on the first pair of timeseries
   [static12, dynamic12] = hoco_fbglmfit(X_full, Obs_full, WSize, param11, param12, param21, TR);

   % HOCo on the second pair of timeseries
   [static13, dynamic13] = hoco_fbglmfit(X_full, T_full, WSize, param11, param12, param21, TR);

   % HOCo on the third pair of timeseries
   [static23, dynamic23] = hoco_fbglmfit(Obs_full, T_full, WSize, param11, param12, param21, TR);

   % ========================================
   % before running kmeans, trim the HOCo results, remove the FC matrices obtained from
   % shrinking/growing windows , because the truncated HOCo results were
   % better than the HOCo itself
   % DON'T CUT THE SHRINKING/GROWING PARTS OF HOCO, ORIGINAL HOCO IS MORE
   % ACCURATE
%    dynamic12.bb  = dynamic12.bb(WSize/2 : size(X_full) - WSize/2);
%    dynamic13.bb  = dynamic13.bb(WSize/2 : size(X_full) - WSize/2);
%    dynamic23.bb  = dynamic23.bb(WSize/2 : size(X_full) - WSize/2);
% %    ========================================
   save(fullfile(Basefile, 'HOCo_results'), 'dynamic12', 'dynamic13', 'dynamic23')
   % kmeans
   kmeans_data_HOCo = [nonzeros(dynamic12.bb), nonzeros(dynamic13.bb), nonzeros(dynamic23.bb)];
   [StateIdx_HOCo, Centroid_HOCo, SumDistPartial_HOCo] = kmeans(kmeans_data_HOCo, k, 'Distance', 'cityblock','MaxIter',100, 'Replicates', 100);  %k=3, for three clusters/states
   % record the time of transition between states
   StateFlip_HOCo{iter} = find(logical(diff(StateIdx_HOCo)));
   % -------------
   % ========================================

end

% Save results to a .mat file
% save('correlated_sinusoids_results.mat', 'all_x', 'all_y', 'all_z', 'all_corr_over_time_xy', 'all_corr_over_time_xz', 'all_corr_over_time_yz');


