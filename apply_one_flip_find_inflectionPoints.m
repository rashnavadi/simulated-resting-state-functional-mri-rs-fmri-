%% August 2024, By Tahereh Rashnavadi
% combine the simulation and finding inflection points in each iterations
clear;
% close all;
Basefile= '/Users/trashnavadi/Documents/Data_Analysis/2022/analyses/kmeans_investigation/2024/august';

% Initialize variables
% Number of iterations
numIterations = 10;
nSubj = numIterations;
% Timepoints
nTimePts = 150;
nStates = 2;
TR = 2; % in sec
time = linspace(0, (nTimePts-1) * TR, nTimePts);

% in seconds
% transitionPoint = [15, 45, 75, 120]; % in TR
transitionPoint = 15;  % Transition between states happens at this point in TR


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

% Initialize a cell array to store the results for each iteration
all_max_values_HOCo = cell(numIterations, 1);
all_max_values_SWC = cell(numIterations, 1);

% Initialize storage for noisy time series
all_noisy_x = zeros(nTimePts, numIterations);
all_noisy_y = zeros(nTimePts, numIterations);
all_noisy_z = zeros(nTimePts, numIterations);


% Run the function for 1000 iterations
for iter = 1:numIterations
%     [x, y, z, corr_over_time_xy, corr_over_time_xz, corr_over_time_yz] = generateCorrelatedSinusoids(nTimePts, transitionPoint);
    [x, y, z] = generateCorrelatedSinusoids(nTimePts, transitionPoint);


    % generate pink noise independently to each timeseries,
    % Generate pink noise for x, y, and z
    pink_noise_x = pinknoise(nTimePts, 1);
    pink_noise_y = pinknoise(nTimePts, 1);
    pink_noise_z = pinknoise(nTimePts, 1);


    % i should keep the ratio between the min and max of the timeseries related
    % to the pink noise the same before and after rescaling
    % Current min and max of the sine signal
    x_noise_min_ratio = min(x)/min(pink_noise_x);
    x_noise_max_ratio = max(x)/max(pink_noise_x);
    x_min_scale = x_noise_min_ratio * min(desired_std/std(pink_noise_x) * pink_noise_x)/ min(x);
    x_max_scale = x_noise_max_ratio * max(desired_std/std(pink_noise_x) * pink_noise_x)/ max(x);
    x_scale = mean([x_min_scale, x_max_scale])/8;

    y_noise_min_ratio = min(y)/min(pink_noise_y);
    y_noise_max_ratio = max(y)/max(pink_noise_y);
    y_min_scale = y_noise_min_ratio * min(desired_std/std(pink_noise_y) * pink_noise_y)/ min(y);
    y_max_scale = y_noise_max_ratio * max(desired_std/std(pink_noise_y) * pink_noise_y)/ max(y);
    y_scale = mean([y_min_scale, y_max_scale])/8;

    z_noise_min_ratio = min(z)/min(pink_noise_z);
    z_noise_max_ratio = max(z)/max(pink_noise_z);
    z_min_scale = z_noise_min_ratio * min(desired_std/std(pink_noise_z) * pink_noise_z)/ min(z);
    z_max_scale = z_noise_max_ratio * max(desired_std/std(pink_noise_z) * pink_noise_z)/ max(z);
    z_scale = mean([z_min_scale, z_max_scale])/8;


    noisy_x = 10000 + (x_scale * x) + (desired_std/std(pink_noise_x) * pink_noise_x);
    noisy_y = 10000 + (y_scale * y) + (desired_std/std(pink_noise_y) * pink_noise_y);
    noisy_z = 10000 + (z_scale * z) + (desired_std/std(pink_noise_z) * pink_noise_z);

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

    %% we are finding the inflections to find the maximum difference in

    % HOCo correlations timeseries of xy, xz and yz
    HOCo_betas = [dynamic12.bb, dynamic13.bb, dynamic23.bb];
    
    % for SWC, r values, correlations timeseries of xy, xz and yz
    SWC_r_values = [FC_21, FC_31, FC_32];

    % Step 1: Compute the derivative of the correlation plots
    delta_HOCo_betas = diff(HOCo_betas);
    delta_SWC_r_values = diff(SWC_r_values);

    %     figure; plot(delta_HOCo_betas(:,1)); hold on; plot(delta_HOCo_betas(:,2)); plot(delta_HOCo_betas(:,3))
    %     title('HOCo derivative of correlations')
    %     figure; plot(delta_SWC_r_values(:,1)); hold on; plot(delta_SWC_r_values(:,2)); plot(delta_SWC_r_values(:,3))
    %     title('SWC derivative of correlations')
        
    % for HOCo, beta values
        % Get the number of detected flips for the current iteration
    num_flips_hoco = numel(StateFlip_HOCo{iter, 1});
    
    % Preallocate an array to store max values for the current iteration
    max_values_HOCo_for_iter = zeros(num_flips_hoco, 1);

    for hoco_detected_flips = 1:numel(StateFlip_HOCo{iter, 1})
        transitionpoint_HOCo = StateFlip_HOCo{iter, 1}(hoco_detected_flips);

        % Step 2: Find the forward inflection point
        % HOCo specific variables
        forward_inflection_1_HOCo = transitionpoint_HOCo;
        xy_inflection_HOCo = delta_HOCo_betas(:, 1);
        for jh = transitionpoint_HOCo:size(xy_inflection_HOCo)
            if jh < size(xy_inflection_HOCo, 1) && (xy_inflection_HOCo(jh) * xy_inflection_HOCo(jh + 1) < 0 || xy_inflection_HOCo(jh) == 0)
                forward_inflection_1_HOCo = jh + 1;
                break;
            end
        end

        forward_inflection_2_HOCo = transitionpoint_HOCo;
        xz_inflection_HOCo = delta_HOCo_betas(:, 2);
        for kh = transitionpoint_HOCo:size(delta_HOCo_betas, 1)
            if kh < size(xz_inflection_HOCo, 1) && (xz_inflection_HOCo(kh) * xz_inflection_HOCo(kh + 1) < 0 || xz_inflection_HOCo(kh) == 0)
                forward_inflection_2_HOCo = kh + 1;
                break;
            end
        end

        forward_inflection_3_HOCo = transitionpoint_HOCo;
        yz_inflection_HOCo = delta_HOCo_betas(:, 3);
        for lh = transitionpoint_HOCo:size(delta_HOCo_betas, 1)
            if lh < size(yz_inflection_HOCo, 1) && (yz_inflection_HOCo(lh) * yz_inflection_HOCo(lh + 1) < 0 || yz_inflection_HOCo(lh) == 0)
                forward_inflection_3_HOCo = lh + 1;
                break;
            end
        end

        % Step 3: Find the backward inflection point for HOCo
        backward_inflection_1_HOCo = transitionpoint_HOCo;
        for jh = transitionpoint_HOCo-1:-1:2
            if (xy_inflection_HOCo(jh - 1) * xy_inflection_HOCo(jh) < 0 || xy_inflection_HOCo(jh - 1) == 0)
                backward_inflection_1_HOCo = jh;
                break;
            end
        end

        backward_inflection_2_HOCo = transitionpoint_HOCo;
        for kh = transitionpoint_HOCo-1:-1:2
            if (xz_inflection_HOCo(kh - 1) * xz_inflection_HOCo(kh) < 0 || xz_inflection_HOCo(kh - 1) == 0)
                backward_inflection_2_HOCo = kh;
                break;
            end
        end

        backward_inflection_3_HOCo = transitionpoint_HOCo;
        for lh = transitionpoint_HOCo-1:-1:2
            if (yz_inflection_HOCo(lh - 1) * yz_inflection_HOCo(lh) < 0 || yz_inflection_HOCo(lh - 1) == 0)
                backward_inflection_3_HOCo = lh;
                break;
            end
        end


        % Step 4: Calculate the correlation values HOCo at the inflection points
        correlation_forward_1_HOCo = HOCo_betas(forward_inflection_1_HOCo, 1);
        correlation_backward_1_HOCo = HOCo_betas(backward_inflection_1_HOCo, 1);

        correlation_forward_2_HOCo = HOCo_betas(forward_inflection_2_HOCo, 2);
        correlation_backward_2_HOCo = HOCo_betas(backward_inflection_2_HOCo, 2);

        correlation_forward_3_HOCo = HOCo_betas(forward_inflection_3_HOCo, 3);
        correlation_backward_3_HOCo = HOCo_betas(backward_inflection_3_HOCo, 3);

        % Step 5: Calculate the correlation values at the inflection points
        % calculate the difference of correlation between the forward and backward inflection points around the
        % transitionpoint for each of the xy, xz and yz timeseries
        % for xy_inflection
        xy_corr_diff_HOCo = correlation_forward_1_HOCo - correlation_backward_1_HOCo;
        xz_corr_diff_HOCo = correlation_forward_2_HOCo - correlation_backward_2_HOCo;
        yz_corr_diff_HOCo = correlation_forward_3_HOCo - correlation_backward_3_HOCo;


        % Step 6:find the maximum of difference between the correlations of the two inflection points around the
        % transitionpoints
        abs_xy_corr_diff_HOCo = abs(xy_corr_diff_HOCo);
        abs_xz_corr_diff_HOCo = abs(xz_corr_diff_HOCo);
        abs_yz_corr_diff_HOCo = abs(yz_corr_diff_HOCo);

        % Find the maximum correlation difference
        max_xy_corr_diff_HOCo = max(abs_xy_corr_diff_HOCo(:));
        max_xz_corr_diff_HOCo = max(abs_xz_corr_diff_HOCo(:));
        max_yz_corr_diff_HOCo = max(abs_yz_corr_diff_HOCo(:));

        % Extracting the values
        val_max_xy_corr_diff_HOCo = full(max_xy_corr_diff_HOCo(1, 1));
        val_max_xz_corr_diff_HOCo = full(max_xz_corr_diff_HOCo(1, 1));
        val_max_yz_corr_diff_HOCo = full(max_yz_corr_diff_HOCo(1, 1));

        % Finding the maximum value among them
        max_value_HOCo = max([val_max_xy_corr_diff_HOCo, val_max_xz_corr_diff_HOCo, val_max_yz_corr_diff_HOCo]);

        % Store the max value for the current flip
        max_values_HOCo_for_flip(hoco_detected_flips) = max_value_HOCo;
    end
    % Store the array of max values for the current iteration
    all_max_values_HOCo{iter} = max_values_HOCo_for_flip;

    % Save the results for this iteration
%     all_max_values_HOCo{iter} = max_values_HOCo;

    % SWC. 
    % Get the number of detected flips for the current iteration
    num_flips_swc = numel(StateFlip_SWC{iter, 1});
   
    % Preallocate an array to store max values for the current iteration
    max_values_SWC_for_flip = zeros(num_flips_swc, 1);

    for swc_detected_flips = 1:numel(StateFlip_SWC{iter, 1})
        transitionpoint_SWC = StateFlip_SWC{iter, 1}(swc_detected_flips);

        % Step 2: Find the forward inflection point
        % SWC specific variables
        forward_inflection_1_SWC = transitionpoint_SWC;
        xy_inflection_SWC = delta_SWC_r_values(:, 1);
        for js = transitionpoint_SWC:size(xy_inflection_SWC)
            if js < size(xy_inflection_SWC, 1) && (xy_inflection_SWC(js) * xy_inflection_SWC(js + 1) < 0 || xy_inflection_SWC(js) == 0)
                forward_inflection_1_SWC = js + 1;
                break;
            end
        end

        forward_inflection_2_SWC = transitionpoint_SWC;
        xz_inflection_SWC = delta_SWC_r_values(:, 2);
        for ks = transitionpoint_SWC:size(delta_SWC_r_values, 1)
            if ks < size(xz_inflection_SWC, 1) && (xz_inflection_SWC(ks) * xz_inflection_SWC(ks + 1) < 0 || xz_inflection_SWC(ks) == 0)
                forward_inflection_2_SWC = ks + 1;
                break;
            end
        end

        forward_inflection_3_SWC = transitionpoint_SWC;
        yz_inflection_SWC = delta_SWC_r_values(:, 3);
        for ls = transitionpoint_SWC:size(delta_SWC_r_values, 1)
            if ls < size(yz_inflection_SWC, 1) && (yz_inflection_SWC(ls) * yz_inflection_SWC(ls + 1) < 0 || yz_inflection_SWC(ls) == 0)
                forward_inflection_3_SWC = ls + 1;
                break;
            end
        end

        % Step 3: Find the backward inflection point for HOCo
        backward_inflection_1_SWC = transitionpoint_SWC;
        for js = transitionpoint_SWC-1:-1:2
            if (xy_inflection_SWC(js - 1) * xy_inflection_SWC(js) < 0 || xy_inflection_SWC(js - 1) == 0)
                backward_inflection_1_SWC = js;
                break;
            end
        end

        backward_inflection_2_SWC = transitionpoint_SWC;
        for ks = transitionpoint_SWC-1:-1:2
            if (xz_inflection_SWC(ks - 1) * xz_inflection_SWC(ks) < 0 || xz_inflection_SWC(ks - 1) == 0)
                backward_inflection_2_SWC = ks;
                break;
            end
        end

        backward_inflection_3_SWC = transitionpoint_SWC;
        for ls = transitionpoint_SWC-1:-1:2
            if (yz_inflection_SWC(ls - 1) * yz_inflection_SWC(ls) < 0 || yz_inflection_SWC(ls - 1) == 0)
                backward_inflection_3_SWC = ls;
                break;
            end
        end

        % Step 4: Calculate the correlation values SWC at the inflection points
        correlation_forward_1_SWC = SWC_r_values(forward_inflection_1_SWC, 1);
        correlation_backward_1_SWC = SWC_r_values(backward_inflection_1_SWC, 1);

        correlation_forward_2_SWC = SWC_r_values(forward_inflection_2_SWC, 2);
        correlation_backward_2_SWC = SWC_r_values(backward_inflection_2_SWC, 2);

        correlation_forward_3_SWC = SWC_r_values(forward_inflection_3_SWC, 3);
        correlation_backward_3_SWC = SWC_r_values(backward_inflection_3_SWC, 3);

        % Step 5: Calculate the correlation values at the inflection points
        % calculate the difference of correlation between the forward and backward inflection points around the
        % transitionpoint for each of the xy, xz and yz timeseries
        % for xy_inflection
        xy_corr_diff_SWC = correlation_forward_1_SWC - correlation_backward_1_SWC;
        xz_corr_diff_SWC = correlation_forward_2_SWC - correlation_backward_2_SWC;
        yz_corr_diff_SWC = correlation_forward_3_SWC - correlation_backward_3_SWC;

        % Step 6: find the maximum of difference between the correlations of the two inflection points around the
        % transitionpoints
        abs_xy_corr_diff_SWC = abs(xy_corr_diff_SWC);
        abs_xz_corr_diff_SWC = abs(xz_corr_diff_SWC);
        abs_yz_corr_diff_SWC = abs(yz_corr_diff_SWC);

        % Find the maximum correlation difference
        max_xy_corr_diff_SWC = max(abs_xy_corr_diff_SWC(:));
        max_xz_corr_diff_SWC = max(abs_xz_corr_diff_SWC(:));
        max_yz_corr_diff_SWC = max(abs_yz_corr_diff_SWC(:));

        % Extracting the values
        val_max_xy_corr_diff_SWC = full(max_xy_corr_diff_SWC(1, 1));
        val_max_xz_corr_diff_SWC = full(max_xz_corr_diff_SWC(1, 1));
        val_max_yz_corr_diff_SWC = full(max_yz_corr_diff_SWC(1, 1));

        % Finding the maximum value among them
        max_value_SWC = max([val_max_xy_corr_diff_SWC, val_max_xz_corr_diff_SWC, val_max_yz_corr_diff_SWC]);

        % Store the max value for the current flip
        max_values_SWC_for_flip(swc_detected_flips) = max_value_SWC;
    end
    % Store the array of max values for the current iteration
    all_max_values_SWC{iter} = max_values_SWC_for_flip;
end
% Display all the HOCo results
disp('HOCo - The maximum values for each detected flip for each iteration are:');
disp(all_max_values_HOCo);

% Display all the SWC results
disp('SWC - The maximum values for each detected flip for each iteration are:');
disp(all_max_values_SWC);

% Save results to a .mat file
% save('one_flip_at_120TR_WS40TR.mat', 'all_noisy_x', 'all_noisy_y', 'all_noisy_z', 'StateFlip_HOCo', 'StateFlip_SWC', 'max_value_HOCo', 'max_value_SWC');



