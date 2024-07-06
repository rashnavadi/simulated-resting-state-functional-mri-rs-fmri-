%% July 2024, By Tahereh Rashnavadi
clear;
close all;

% Initialize variables
% Number of iterations
numIterations = 10;
nSubj = numIterations;
% Timepoints
nTimePts = 150;
nStates = 2;
TR = 2; % in sec
transitionPoint = 75;  % Transition between states happens at this point
k = nStates;
% WSize = [38, 50, 76, 100];
WSize = 50; % in timepoints (TRs), 300 sec
window_size = WSize;
% Window size for sliding window correlation
window_size = WSize;
half_window = floor(window_size / 2);
% generate three timeseries of x, y and z with desired correlations between
% them: corr(x, y(1:transitionPoint))=.3, corr(x, y(transitionPoint+1:end)=.8; corr(y, z) =.8, corr(x, z(1:transitionPoint)=.8,
% corr(x, z(transitionPoint+1:end)=.3
% transitionPoint=[15, 45, 75, 120]; in TR 
transitionPoint = 120;
t = 0:2:298;
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

% you can generate 1/f noise (Called pinknoise) in MATLAB directly.
% For a sampling rate of 0.5 HZ (TR=2sec)and a duration of 300 seconds (150 TR):
fs = 0.5;
duration = 300;
N = fs * duration; % number of samples
pink_noise_series = cell(3, numIterations);
% pinknoise(0.5 * 300); % in the for loop to generate different sets of
% noise in each iteration

% We should also first add Gaussian noise to the timeseries itself if the high-frequency pinknoise goes to zero:
% Where I is the timeseries, m is the mean, and var_gauss is the
% variance. J = imnoise(I,'gaussian',m,var_gauss)
% You want to end up with a signal that has a signal-to-noise ratio of about 50, ...
% which would approximate a real scenario of an average over an ROI.
% default valeus

% Preallocate arrays to store state flips
StateFlip_SWC = cell(numIterations, 1);

% Initialize storage for noisy time series
noisy_x = zeros(length(x), numIterations);
noisy_y = zeros(length(y), numIterations);
noisy_z = zeros(length(z), numIterations);
all_x = zeros(original_length, numIterations);
all_y = zeros(original_length, numIterations);
all_z = zeros(original_length, numIterations);

% Run the function for 1000 iterations
for iter = 1:numIterations
    % generate pink noise independently to each timeseries, otherwise
    % for large noise levels rather than making the TSs more
    % uncorrelated they become more alike and highly correlated
    pink_noise_x = pinknoise(N);
    pink_noise_y = pinknoise(N);
    pink_noise_z = pinknoise(N);
    % Store the pink noise in the cell array
    pink_noise_series{iter} = [pink_noise_x, pink_noise_y, pink_noise_z] ;

    % add gaussian and pink noise
    noisy_x = imnoise(x, 'gaussian', 0, .01) + pink_noise_x;
    noisy_y = imnoise(y, 'gaussian', 0, .01) + pink_noise_y;
    noisy_z = imnoise(z, 'gaussian', 0, .01) + pink_noise_z;
    % while adding the Gaussian noise, desired SNR in dB (SNR of each timeseries after adding noise should be that)
    % The SNR is 50, not 50 dB: In fMRI lingo, SNR = mean signal over time / variance of noise over time
    %Calculate the spatial SNR: spatial_SNR = mean_signal / noise_std;
    SNR_fmri = 50;
    % then the std of the noise to meet this snr, Desired mean and standard deviation
    desired_mean = 10000;
%     desired_std = 200;
    % Adjust noisy_x
    current_mean_x = mean(noisy_x) + 10000; % signal mean scaled by adding 10000
    current_std_noise_x = std(noisy_x - x); % noise std
    snr_noisy_x = current_mean_x/current_std_noise_x;
    % Adjust noisy_y
    current_mean_y = mean(noisy_y) + 10000;
    current_std_noise_y = std(noisy_y - y);
    snr_noisy_y = current_mean_y/current_std_noise_y;
    % Adjust noisy_z
    current_mean_z = mean(noisy_z) + 10000;
    current_std_noise_z = std(noisy_z - z);
    snr_noisy_z = current_mean_z/current_std_noise_z;

    % Check condition and store the noisy_x if condition is met
    if current_mean_x/current_std_noise_x >= 50 && current_mean_y/current_std_noise_y >= 50 && current_mean_z/current_std_noise_z >= 50
        noisy_x_set = noisy_x;
        noisy_y_set = noisy_y;
        noisy_z_set = noisy_z;
    end

    % Store the result
    all_x(:, iter) = noisy_x_set;
    all_y(:, iter) = noisy_y_set;
    all_z(:, iter) = noisy_z_set;

    % Sliding window cross-correlation
    num_windows = original_length - window_size + 1;
    swc_matrix = zeros(num_windows, 3, 3); % Assuming 3 time series (x, y, z)
    
    for w = 1:num_windows
        window_x = noisy_x_set(w:w + window_size - 1);
        window_y = noisy_y_set(w:w + window_size - 1);
        window_z = noisy_z_set(w:w + window_size - 1);
        
        corr_matrix = corr([window_x, window_y, window_z]);
        swc_matrix(w, :, :) = corr_matrix;
    end
    
    % Flatten the sliding window correlation matrices for clustering
    swc_flattened = reshape(swc_matrix, [num_windows, 9]);
    
    % K-means clustering
    [idx, C] = kmeans(swc_flattened, k);
    
    % Store state flips
    StateFlip_SWC{iter} = idx;
end

% Save the results to a .mat file
save('timeseries_with_noise.mat', 'all_x', 'all_y', 'all_z', 'StateFlip_SWC');