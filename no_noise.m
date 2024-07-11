%% July 2024, By Tahereh Rashnavadi
clear;
% close all;
Basefile= '/Users/trashnavadi/Documents/Data_Analysis/2022/analyses/kmeans_investigation/2024/July';

% Initialize variables
% Number of iterations
numIterations = 1000;
nSubj = numIterations;
% Timepoints
nTimePts = 150;
nStates = 2;
TR = 2; % in sec
transitionPoint = 75;  % Transition between states happens at this point
k = nStates;
% WSize = [38, 50, 76, 100];
WSize = 50; % in timepoints (TRs), 300 sec

% Load the original noise-free time series from a file (B-G)
% impose two states on x, y and z timeseries
% Generate time vector
t = 0:TR:(TR * nTimePts - 2);
fs = 0.06; % frequency for resting state fmri
x = sin(2 * pi * t * 0.06);
x = x';
x = zscore(x);
y_low  = sin(2 * pi * (t(1:transitionPoint) * 0.06 + 0.1));  
y_high = sin(2 * pi * (t(transitionPoint+1:end) * 0.06 + 0.2));
y = [y_low, y_high];
y = y';
y = zscore(y);
z_low  = sin(2 * pi * (t(1:transitionPoint) * 0.06 + 0.2));
z_high = sin(2 * pi * (t(transitionPoint+1:end) * 0.06 + 0.1));
z = [z_low, z_high];
z = z';
z = zscore(z);
original_length = nTimePts;

% Preallocate arrays to store state flips
StateFlip_SWC = cell(numIterations, 1);
StateFlip_HOCo = cell(numIterations, 1);

% Initialize storage for 1000 samples sets of xyz
all_x = zeros(original_length, numIterations);
all_y = zeros(original_length, numIterations);
all_z = zeros(original_length, numIterations);

% Run the function for 1000 iterations
for iter = 1:numIterations
    % Add random noise to the time series
    noise_level = 0.001; % Adjust the noise level as needed
    noisy_x_set = x + noise_level * randn(size(x));
    noisy_y_set = y + noise_level * randn(size(y));
    noisy_z_set = z + noise_level * randn(size(z));

    % Store the result
    all_x(:, iter) = noisy_x_set;
    all_y(:, iter) = noisy_y_set;
    all_z(:, iter) = noisy_z_set;

    % Prepare data for HOCo and SWC
    TSs = cell(2, 3);
    TSs{1, 1} = 'timeseries_1';
    TSs{2, 1} = all_x(:, iter);
    TSs{1, 2} = 'timeseries_2';
    TSs{2, 2} = all_y(:, iter);
    TSs{1, 3} = 'timeseries_3';
    TSs{2, 3} = all_z(:, iter);
    save(fullfile(Basefile, 'sim_timeseries_1.mat'), 'TSs')

    % HOCo and SWC processing
    load(fullfile(Basefile, 'sim_timeseries_1.mat'))
    SaveSuffix = 'Pearsons'; SurrOption = 0;
    Time_series_file = fullfile(Basefile, 'sim_timeseries_1.mat');
    hoco_FunConn(Time_series_file, WSize, SurrOption, SaveSuffix)
    load(fullfile(Basefile, 'mtx_SldWFCPearsons.mat'))
    FC_21 = SldWFCArray{3, 2};
    FC_31 = SldWFCArray{4, 2};
    FC_32 = SldWFCArray{4, 3};

    % SWC Analysis
    kmeans_data_SWC = [FC_21, FC_31, FC_32];
    [StateIdx_SWC, Centroid_SWC, SumDistPartial_SWC] = kmeans(kmeans_data_SWC, k, 'Distance', 'cityblock', 'MaxIter', 100, 'Replicates', 100);
    StateFlip_SWC{iter, 1} = find(logical(diff(StateIdx_SWC)));

    % HOCo Analysis
    X_full = TSs{2, 1};
    Obs_full = TSs{2, 2};
    T_full = TSs{2, 3};
    param11 = 1; param12 = 1; param21 = 1;

    [static12, dynamic12] = hoco_fbglmfit(X_full, Obs_full, WSize, param11, param12, param21, TR);
    [static13, dynamic13] = hoco_fbglmfit(X_full, T_full, WSize, param11, param12, param21, TR);
    [static23, dynamic23] = hoco_fbglmfit(Obs_full, T_full, WSize, param11, param12, param21, TR);

    % ========================================
    %     % before running kmeans, trim the HOCo results, remove the FC matrices obtained from
    % shrinking/growing windows , because the truncated HOCo results were
    % better than the HOCo itself
%         dynamic12.bb  = dynamic12.bb(WSize/2 : size(X_full) - WSize/2);
%         dynamic13.bb  = dynamic13.bb(WSize/2 : size(X_full) - WSize/2);
%         dynamic23.bb  = dynamic23.bb(WSize/2 : size(X_full) - WSize/2);

    kmeans_data_HOCo = [nonzeros(dynamic12.bb), nonzeros(dynamic13.bb), nonzeros(dynamic23.bb)];
    [StateIdx_HOCo, Centroid_HOCo, SumDistPartial_HOCo] = kmeans(kmeans_data_HOCo, k, 'Distance', 'cityblock', 'MaxIter', 100, 'Replicates', 100);
    StateFlip_HOCo{iter, 1} = find(logical(diff(StateIdx_HOCo)));
end

% Save the results to a .mat file
save('timeseries_with_noise.mat', 'all_x', 'all_y', 'all_z', 'StateFlip_HOCo', 'StateFlip_SWC');




