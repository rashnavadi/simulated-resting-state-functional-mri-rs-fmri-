function [noisy_x_set, noisy_y_set, noisy_z_set, corr_xy, corr_xz, corr_yz, corr_over_time_xy, corr_over_time_xz, corr_over_time_yz ] = generateCorrelatedSinusoids(nTimePts, transitionPoint)
    
    TR = 2; % Time resolution in seconds

    ntimepoints = nTimePts;

    % Generate time vector
    t = 0:TR:(TR * ntimepoints - 2);

    % Generate the first timeseries (x: a sine wave with period t/2)
    x = sin(2 * pi * t * 0.06);

    % Desired correlations for y
    desired_corr_xy_low = 0.3;
    desired_corr_xy_high = 0.8;

    % Desired correlations for z
    desired_corr_xz = 0.8;
    desired_corr_yz_low = 0.8;
    desired_corr_yz_high = 0.3;

    % Adjust the phase of the second sine wave to achieve the desired correlation
    phi_xy_low = acos(desired_corr_xy_low);  % Phase shift to achieve desired correlation
    phi_xy_high = acos(desired_corr_xy_high);  % Phase shift to achieve desired correlation

    % Generate the second timeseries (y, another sine wave)
    y = zeros(size(x));
    y(1:transitionPoint) = sin(2 * pi * (t(1:transitionPoint) * 0.06) + phi_xy_low);
    y(transitionPoint+1:end) = sin(2 * pi * (t(transitionPoint+1:end) * 0.06) + phi_xy_high);

    % Generate the third timeseries (z, another sine wave with changing correlation to y)
    z = zeros(size(x));
    z(1:transitionPoint) = sin(2 * pi * (t(1:transitionPoint) * 0.06) + acos(desired_corr_xz));
    z(transitionPoint+1:end) = sin(2 * pi * (t(transitionPoint+1:end) * 0.06) + acos(desired_corr_yz_high));

    % Add random noise to the time series
    noise_level = 0.001 ; % Adjust the noise level as needed
    noisy_x_set = x + noise_level * randn(size(x));
    noisy_y_set = y + noise_level * randn(size(y));
    noisy_z_set = z + noise_level * randn(size(z));
    
    % Verify correlations
    corr_xy = corr(noisy_x_set', noisy_y_set');
    corr_xz = corr(noisy_x_set', noisy_z_set');
    corr_yz = corr(noisy_y_set', noisy_z_set');

    disp(['Correlation between x and z: ', num2str(corr_xz)]);
    disp(['Correlation between y and z: ', num2str(corr_yz)]);

    % Plot the timeseries to visualize the abrupt change
%     figure;
%     subplot(3,1,1);
%     plot(t, x, 'b', 'DisplayName', 'x');
%     hold on;
%     plot(t, y, 'r', 'DisplayName', 'y');
%     plot(t, z, 'g', 'DisplayName', 'z');
%     xlabel('Time (s)');
%     ylabel('Amplitude');
%     title('Time Series');
%     legend;
%     ylim([-0.01, 1.01]); % Set y-axis limits


%     % Calculate the correlation over time
    window_size = 50; % Size of the sliding window for correlation calculation: 30, 40, 50
%     corr_over_time_xy = zeros(1, nTimePts - window_size + 1);
%     corr_over_time_xz = zeros(1, nTimePts - window_size + 1);
%     corr_over_time_yz = zeros(1, nTimePts - window_size + 1);
% 
    for i = 1:nTimePts - window_size + 1
        corr_over_time_xy(i) = corr(noisy_x_set(i:i+window_size-1)', noisy_y_set(i:i+window_size-1)');
        corr_over_time_xz(i) = corr(noisy_x_set(i:i+window_size-1)', noisy_z_set(i:i+window_size-1)');
        corr_over_time_yz(i) = corr(noisy_y_set(i:i+window_size-1)', noisy_z_set(i:i+window_size-1)');
    end
% 
%     subplot(3,1,2);
%     plot(t(window_size:end), corr_over_time_xy, 'k', 'DisplayName', 'Correlation xy');
%     xlabel('Time (s)');
%     ylabel('Correlation');
%     title('Correlation between x and y over Time');
%     legend;
%     ylim([-0.01, 1.01]); % Set y-axis limits
% 
%     subplot(3,1,3);
%     plot(t(window_size:end), corr_over_time_xz, 'm', 'DisplayName', 'Correlation xz');
%     hold on;
%     plot(t(window_size:end), corr_over_time_yz, 'c', 'DisplayName', 'Correlation yz');
%     xlabel('Time (s)');
%     ylabel('Correlation');
%     title('Correlation between x and z and y and z over Time');
%     legend;
%     ylim([-0.01, 1.01]); % Set y-axis limits

    % Transpose the outputs to match the expected output format
    noisy_x_set = noisy_x_set';
    noisy_y_set = noisy_y_set';
    noisy_z_set = noisy_z_set';

end
