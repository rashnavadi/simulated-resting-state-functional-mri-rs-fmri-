%% July 2024, By Tahereh Rashnavadi
clear;
close all;
Basefile= '/Users/trashnavadi/Documents/Data_Analysis/2022/analyses/kmeans_investigation/2024/July';

% Initialize variables
% Number of iterations
numIterations = 1000;
nSubj = numIterations;
% Timepoints
nTimePts = 150;
nStates = 2;
TR = 2; % in sec
time = linspace(0, (nTimePts-1) * TR, nTimePts);

% in seconds
% transitionPoint = [15, 45, 75, 120]; % in TR
transitionPoint = 15;  % Transition between states happens at this point in TR
% transitionPoint1 = 15;  % Transition between states happens at this point
% transitionPoint2 = 15;  % Transition between states happens at this point

k = nStates;
% WSize = [30, 40, 50]; % in timepoints (TR = 2sec)
WSize = 50; % in timepoints (TR = 2sec)

% you can generate 1/f noise (Called pinknoise) in MATLAB directly.
% For a sampling rate of 0.5 HZ (TR=2sec)and a duration of 300 seconds (150 TR):
fs = 0.5;
duration = 300;
% generate pinknoise
pink_noise_obj = dsp.ColoredNoise('Color', 'pink', 'SamplesPerFrame', nTimePts, 'NumChannels', 1);
pink_noise_series = cell(3, numIterations);
desired_snr = 50;
desired_std = 200;


% Preallocate arrays to store state flips
StateFlip_SWC = cell(numIterations, 1);
StateFlip_HOCo = cell(numIterations, 1);

% Initialize storage for noisy time series
all_noisy_x = zeros(nTimePts, numIterations);
all_noisy_y = zeros(nTimePts, numIterations);
all_noisy_z = zeros(nTimePts, numIterations);


% Run the function for 1000 iterations
for iter = 1:numIterations
    [x, y, z, corr_over_time_xy, corr_over_time_xz, corr_over_time_yz] = generateCorrelatedSinusoids(nTimePts, transitionPoint);
    %    [x, y, z, corr_over_time_xy, corr_over_time_xz, corr_over_time_yz] = generateCorrelatedSinusoids_2flips(nTimePts, transitionPoint1, transitionPoint);

    % generate pink noise independently to each timeseries,
    % Generate pink noise for x, y, and z
    pink_noise_x = pinknoise(nTimePts, 1);
    pink_noise_y = pinknoise(nTimePts, 1);
    pink_noise_z = pinknoise(nTimePts, 1);


    % i should keep the ratio between the min and max of the timeseries related
    % to the pink noise the same before and after rescaling
    % Current min and max of the sine signal
    scaling_factor_x_min = min(x) / min(pink_noise_x) * abs(min(desired_std/std(pink_noise_x) * pink_noise_x));
    scaling_factor_x_max = max(x) / max(pink_noise_x) * max(desired_std/std(pink_noise_x) * pink_noise_x);
    scaling_factor_x = max([scaling_factor_x_min, scaling_factor_x_max])/10;

    scaling_factor_y_min = min(y) / min(pink_noise_y) * abs(min(desired_std/std(pink_noise_y) * pink_noise_y));
    scaling_factor_y_max = max(y) / max(pink_noise_y) * max(desired_std/std(pink_noise_y) * pink_noise_y);
    scaling_factor_y = max([scaling_factor_y_min, scaling_factor_y_max])/10;

    scaling_factor_z_min = min(z) / min(pink_noise_z) * abs(min(desired_std/std(pink_noise_z) * pink_noise_z));
    scaling_factor_z_max = max(z) / max(pink_noise_z) * max(desired_std/std(pink_noise_z) * pink_noise_z);
    scaling_factor_z = max([scaling_factor_z_min, scaling_factor_z_max])/10;

    noisy_x = 10000 + (scaling_factor_x * x) + (desired_std/std(pink_noise_x) * pink_noise_x);
    noisy_y = 10000 + (scaling_factor_y * y) + (desired_std/std(pink_noise_y) * pink_noise_y);
    noisy_z = 10000 + (scaling_factor_z * z) + (desired_std/std(pink_noise_z) * pink_noise_z);

%     noisy_x = zscore(noisy_x);
%     noisy_y = zscore(noisy_y);
%     noisy_z = zscore(noisy_z);

    % Store the noisy signals
    all_noisy_x(:, iter) = noisy_x;
    all_noisy_y(:, iter) = noisy_y;
    all_noisy_z(:, iter) = noisy_z;

    % Verify the standard deviation of the noise
    final_noise_x = noisy_x - x;
    final_noise_y = noisy_y - y;
    final_noise_z = noisy_z - z;
    final_std_noise_x = std(final_noise_x);
    final_std_noise_y = std(final_noise_y);
    final_std_noise_z = std(final_noise_z);

    % Display the result
    disp(['Final standard deviation of the noise for x: ', num2str(final_std_noise_x)]);
    disp(['Final standard deviation of the noise for y: ', num2str(final_std_noise_y)]);
    disp(['Final standard deviation of the noise for z: ', num2str(final_std_noise_z)]);

    % Prepare data for HOCo and SWC
    TSs = cell(2, 3);
    TSs{1, 1} = 'timeseries_1';
    TSs{2, 1} = all_noisy_x(:, iter);
    TSs{1, 2} = 'timeseries_2';
    TSs{2, 2} = all_noisy_y(:, iter);
    TSs{1, 3} = 'timeseries_3';
    TSs{2, 3} = all_noisy_z(:, iter);
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

% Save results to a .mat file
save('one_flip_at_15TR_WS50TR.mat', 'all_noisy_x', 'all_noisy_y', 'all_noisy_z', 'StateFlip_HOCo', 'StateFlip_SWC');



