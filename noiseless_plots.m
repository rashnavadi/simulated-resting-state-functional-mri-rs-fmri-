%% two states and one flip/transition
% Example data: True Positives for SWC and HOCo
% Data for True Positives (TP) percentages for SWC and HOCo across different acceptance ranges



%% noise-less 2 states 1 transition
% transition at 15 TR
tp_count_SWC = [
    0, 0, 399, 414; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    0, 0, 0, 384;   % for WSize 40
    0, 0, 0, 0      % for WSize 50
];

tp_count_HOCo = [
    0, 0, 0, 453; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    0, 0, 0, 0;   % for WSize 40
    0, 0, 0, 0   % for WSize 50
];

% Mean and std deviation for SWC
mean_SWC = {
    [NaN, NaN, 23.9, 23.98];
    [NaN, NaN, NaN, 28.93];
    [NaN, NaN, NaN, NaN];
};
std_SWC = {
    [NaN, NaN, 0.43, 0.58];
    [NaN, NaN,NaN, 0.49];
    [NaN, NaN, NaN, NaN];
};

% Mean and std deviation for HOCo
mean_HOCo = {
    [NaN, NaN, NaN, 28.72];
    [NaN, NaN, NaN, NaN];
    [NaN, NaN, NaN, NaN]
};
std_HOCo = {
    [NaN, NaN,NaN,0.86];
    [NaN, NaN, NaN, NaN];
    [NaN, NaN, NaN, NaN]
};

% Detection percentage annotations
correct_flip_percentage_SWC = [425, 410, 410]/500*100; % in percentage for three wsizes
correct_flip_percentage_HOCo = [500, 500, 500]/500*100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transition at 45 TR
% True positive counts for each method and window size
tp_count_SWC = [
    [492, 500, 500, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [4, 373, 500, 500];   % for WSize 40
    [0, 0, 372, 500]      % for WSize 50
];

tp_count_HOCo = [
    [475, 500, 500, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [240, 500, 500, 500];   % for WSize 40
    [44, 488, 500, 500]    % for WSize 50
];

% Mean and std deviation for SWC
mean_SWC = {
    [45.65, 45.69, 45.69, 45.69];
    [47, 49.31, 49.77, 49.77];
    [NaN, NaN, 54.42, 54.92]
};
std_SWC = {
    [0.78, 0.83, 0.83, 0.83];
    [0, 0.58, 0.96, 0.96];
    [NaN, NaN, 0.57, 1.04]
};

% Mean and std deviation for HOCo
mean_HOCo = {
    [46.39, 46.47, 46.47, 46.47];
    [46.85, 47.61, 47.61, 47.61];
    [46.95, 48.67, 48.73, 48.73]
};
std_HOCo = {
    [0.68, 0.76, 0.76, 0.76];
    [0.35, 0.85, 0.85, 0.85];
    [0.21, 0.83, 0.91, 0.91]
};

% Detection percentage annotations
correct_flip_percentage_SWC = [100, 100, 100]; % in percentage for three wsizes
correct_flip_percentage_HOCo = [100, 100, 100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transition at 75 TR
tp_count_SWC = [
    [491, 500, 500, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [458, 500, 500, 500];   % for WSize 40
    [450, 500, 500, 500]      % for WSize 50
];

tp_count_HOCo = [
    [485, 500, 500, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [472, 500, 500, 500];   % for WSize 40
    [459, 500, 500, 500]    % for WSize 50
];

% Mean and std deviation for SWC
mean_SWC = {
    [76.29, 76.32, 76.32, 76.32];
    [76.04, 76.22, 76.22, 76.22];
    [75.8, 76.02, 76.02, 76.02]
};
std_SWC = {
    [0.63, 0.67, 0.67, 0.67];
    [0.79, 0.96, 0.96, 0.96];
    [0.79, 1.02, 1.02, 1.02]
};

% Mean and std deviation for HOCo
mean_HOCo = {
    [76.09, 76.15, 76.15, 76.15];
    [76.16, 76.26, 76.26, 76.26];
    [76.09, 76.26, 76.26, 76.26]
};
std_HOCo = {
    [0.7, 0.76, 0.76, 0.76];
    [0.79, 0.88, 0.88, 0.88];
    [0.85, 0.99, 0.99, 0.99]
};

% Detection percentage annotations
correct_flip_percentage_SWC = [500, 500, 500]/500*100; % in percentage for three wsizes
correct_flip_percentage_HOCo = [500, 500, 500]/500*100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transition at 120 TR
tp_count_SWC = [
    [9, 496, 500, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [0, 0, 490, 500];   % for WSize 40
    [0, 0, 0, 475]      % for WSize 50
];

tp_count_HOCo = [
    [298, 500, 500, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [5, 395, 500, 500];   % for WSize 40
    [0, 37, 499, 500]    % for WSize 50
];

% Mean and std deviation for SWC
mean_SWC = {
    [118.11, 116.15, 116.14, 116.14];
    [NaN, NaN, 111.13, 11.09];
    [NaN, NaN, NaN, 106.09]
};
std_SWC = {
    [0.33, 0.94, 0.96, 0.96];
    [NaN, NaN, 0.8, 0.85];
    [NaN, NaN, NaN, 0.94]
};

% Mean and std deviation for HOCo
mean_HOCo = {
    [118.3, 117.73, 117.73, 117.73];
    [118, 115.56, 115.2, 115.2];
    [NaN, 115, 112.85, 112.84]
};
std_HOCo = {
    [0.48, 0.81, 0.81, 0.81];
    [0, 0.69, 0.95, 0.95];
    [NaN, 0, 1.12, 1.13]
};

% Detection percentage annotations
correct_flip_percentage_SWC = [500, 500, 500]/500*100; % in percentage for three wsizes
correct_flip_percentage_HOCo = [500, 500, 500]/500*100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
for i = 1:num_wsize
    subplot(1, num_wsize, i); % Adjust subplots based on number of window sizes
    
    % Smooth line connection for SWC and HOCo using pchip interpolation
    xq = linspace(min(win_ranges), max(win_ranges), 100);
    yq_swc = interp1(win_ranges, tp_count_SWC(i, :), xq, 'pchip'); % pchip interpolation for smooth connection
    yq_hoco = interp1(win_ranges, tp_count_HOCo(i, :), xq, 'pchip');
    
    % Plot the smooth lines for SWC and HOCo
    p1 = plot(xq, yq_swc, '-', 'LineWidth', 1.5, 'Color', 'blue', 'DisplayName', 'SWC');
    hold on;
    p2 = plot(xq, yq_hoco, '-', 'LineWidth', 1.5, 'Color', 'red', 'DisplayName', 'HOCo');
    
    % Plot the original data points with markers and error bars
    for j = 1:length(win_ranges)
        % SWC data points with error bars
        if ~isnan(mean_SWC{i}(j))
            errorbar(win_ranges(j)-0.2, tp_count_SWC(i, j), std_SWC{i}(j), 'o', 'LineWidth', 1.5, ...
                     'MarkerSize', 8, 'Color', 'blue', 'CapSize', 10, 'MarkerFaceColor', 'blue');
        end
        
        % HOCo data points with error bars
        if ~isnan(mean_HOCo{i}(j))
            errorbar(win_ranges(j)+0.2, tp_count_HOCo(i, j), std_HOCo{i}(j), 's', 'LineWidth', 1.5, ...
                     'MarkerSize', 8, 'Color', 'red', 'CapSize', 10, 'MarkerFaceColor', 'red');
        end
    end
    
    % Annotation text for SWC and HOCo correct detection percentages
    swc_annotation = sprintf('SWC Correct Detection: %.1f%%', correct_flip_percentage_SWC(i));
    hoco_annotation = sprintf('HOCo Correct Detection: %.1f%%', correct_flip_percentage_HOCo(i));
    
    % Set x-ticks to only show the acceptance window ranges
    set(gca, 'XTick', win_ranges);
    
    % Only show y-axis label on the first subplot
    if i == 1
        ylabel('True Positives (out of 500)', 'FontSize', 12);
    end
    xlabel('Acceptance Window Range (TR)', 'FontSize', 12);
    
    title(['WSize = ', num2str(wsize(i)), ' TR'], 'FontSize', 14);
    grid off;
    ylim([0, 500]); % Set y-axis limits to 0-500 for clarity
    
    % Add a subplot-specific legend with correct detection percentages for each subplot
    lgd = legend({['SWC (' swc_annotation ')'], ['HOCo (' hoco_annotation ')']}, 'Location', 'southeast');
    lgd.FontSize = 10;
    
    hold off;
end

sgtitle('Transition at 120 TR', 'FontSize', 16);

%% %% noiseless two states and two transitions
% Data for True Positives (TP) percentages for SWC and HOCo for transition imposed at 15 and 45 TR
tp_count_SWC = [
    [0, 0, 0, 0]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [0, 0, 0, 0];   % for WSize 40
    [0, 0, 0, 0]      % for WSize 50
];


tp_count_HOCo = [
    [0, 1, 499, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [0, 0, 0, 0];   % for WSize 40
    [0, 0, 0, 0]   % for WSize 50
];

% Mean and std deviation for SWC
mean1_SWC = {
    [NaN, NaN, NaN, NaN];
    [NaN, NaN, NaN, NaN];
    [NaN, NaN, NaN, NaN];
};
std1_SWC = {
    [NaN, NaN, NaN, NaN];
    [NaN, NaN, NaN, NaN];
    [NaN, NaN, NaN, NaN];
};
mean2_SWC = {
    [48.84, 48.84, 48.84, 48.84];
    [51.41, 51.41, 51.41, 51.41];
    [55, 55, 55, 55, 55];
};
std2_SWC = {
    [0.6, 0.6, 0.6, 0.6];
    [0.57, 0.57, 0.57, 0.57];
    [0.66, 0.66, 0.66, 0.66];
};

% Mean and std deviation for HOCo
mean1_HOCo = {
    [NaN, 10, 7.31, 7.3];
    [NaN, NaN, NaN, NaN];
    [NaN, NaN, NaN, NaN]
};
std1_HOCo = {
    [NaN, 0, 1, 1];
    [NaN, NaN, NaN, NaN];
    [NaN, NaN, NaN, NaN]
};
mean2_HOCo = {
    [NaN, 50, 51.98, 51.98];
    [57.15, 57.15, 57.15, 57.15];
    [60.21, 60.21, 60.21, 60.21]
};
std2_HOCo = {
    [NaN, 0, 0.58, 0.59];
    [0.62, 0.62, 0.62, 0.62];
    [0.72, 0.72, 0.72, 0.72]
};

% Detection percentage annotations
correct_flip_percentage_SWC = [0, 0, 0]/500*100; % in percentage for three wsizes
correct_flip_percentage_HOCo = [500, 0, 0]/500*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data for True Positives (TP) percentages for SWC and HOCo for transition imposed at 15 and 75 TR
% tp_count_SWC = [
%     [0, 0, 0, 0]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
%     [0, 0, 0, 0];   % for WSize 40
%     [0, 0, 0, 0]      % for WSize 50
% ];
% 
% tp_count_HOCo = [
%     [491, 500, 500, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
%     [485, 500, 500, 500];   % for WSize 40
%     [0, 492, 500, 500]   % for WSize 50
% ];
% 
% % Mean and std deviation for SWC
% mean1_SWC = {
%     [NaN, NaN, NaN, NaN];
%     [NaN, NaN, NaN, NaN];
%     [NaN, NaN, NaN, NaN];
% };
% std1_SWC = {
%     [NaN, NaN, NaN, NaN];
%     [NaN, NaN, NaN, NaN];
%     [NaN, NaN, NaN, NaN];
% };
% mean2_SWC = {
%     [76.28, 76.28, 76.28, 76.28];
%     [77.26, 77.26, 77.26, 77.26];
%     [79.59, 79.59, 79.59, 79.59];
% };
% std2_SWC = {
%     [2.74, 2.74, 2.74, 2.74];
%     [0.58, 0.58, 0.58, 0.58];
%     [0.67, 0.67, 0.67, 0.67];
% };
% 
% % Mean and std deviation for HOCo
% mean1_HOCo = {
%     [14.56, 14.55, 14.55, 14.55];
%     [14.9, 14.88, 14.88, 14.88];
%     [NaN, 12.54, 12.5, 12.5]
% };
% std1_HOCo = {
%     [0.65, 0.66, 0.66, 0.66];
%     [0.86, 0.88, 0.88, 0.88 ];
%     [NaN, 1.16, 1.21, 1.21]
% };
% mean2_HOCo = {
%     [76.57, 76.59, 76.59, 76.59];
%     [76.52, 76.56, 76.56, 76.56];
%     [NaN, 79.25, 79.27, 79.27]
% };
% std2_HOCo = {
%     [0.5, 0.53, 0.53, 0.53];
%     [0.52, 0.57, 0.57, 0.57];
%     [NaN, 0.57, 0.59, 0.59]
% };
% 
% % Detection percentage annotations
% correct_flip_percentage_SWC = [1, 0, 0]/500*100; % in percentage for three wsizes
% correct_flip_percentage_HOCo = [500, 500, 500]/500*100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data for True Positives (TP) percentages for SWC and HOCo for transition imposed at 15 and 120 TR
% tp_count_SWC = [
%     [0, 52, 500, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
%     [0, 0, 116, 500];   % for WSize 40
%     [0, 0, 0, 97]      % for WSize 50
% ];
% 
% tp_count_HOCo = [
%     [0, 0, 500, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
%     [0, 0, 261, 500];   % for WSize 40
%     [0, 0, 0, 459]   % for WSize 50
% ];
% 
% % Mean and std deviation for SWC
% mean1_SWC = {
%     [NaN, 24.6, 21.05, 21.05];
%     [NaN, NaN, 24.8, 25.63];
%     [NaN, NaN, NaN, 9.57];
% };
% std1_SWC = {
%     [NaN, 0.5, 0.87, 0.87];
%     [NaN, NaN, 0.4, 0.76];
%     [NaN, NaN, NaN, 0.5];
% };
% mean2_SWC = {
%     [NaN, 120.04, 124.22, 124.22];
%     [NaN, NaN, 110.03, 109.46];
%     [NaN, NaN, NaN, 105.02];
% };
% std2_SWC = {
%     [NaN, 0.19, 0.46, 0.46];
%     [NaN, NaN, 0.18, 0.55];
%     [NaN, NaN, NaN, 0.14];
% };
% 
% % Mean and std deviation for HOCo
% mean1_HOCo = {
%     [NaN, NaN, 21.92, 21.92];
%     [NaN, NaN, 24.93, 25.47];
%     [NaN, NaN, NaN, 29.67]
% };
% std1_HOCo = {
%     [NaN, NaN, 0.46, 0.46];
%     [NaN, NaN, 0.25, 0.6];
%     [NaN, NaN, NaN, 0.48]
% };
% mean2_HOCo = {
%     [NaN, NaN, 114.63, 114.63];
%     [NaN, NaN, 111.15, 111.04];
%     [NaN, NaN, NaN, 106.68]
% };
% std2_HOCo = {
%     [NaN, NaN, 0.56, 0.56];
%     [NaN, NaN, 0.63, 0.62];
%     [NaN, NaN, NaN, 0.63]
% };
% 
% % Detection percentage annotations
% correct_flip_percentage_SWC = [500, 500, 500]/500*100; % in percentage for three wsizes
% correct_flip_percentage_HOCo = [500, 500, 500]/500*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data for True Positives (TP) percentages for SWC and HOCo for transition imposed at 45 and 75 TR
tp_count_SWC = [
    [0, 302, 500, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [0, 0, 499, 500];   % for WSize 40
    [0, 0, 0, 500]      % for WSize 50
];

tp_count_HOCo = [
    [0, 2, 500, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [0, 0, 271, 500];   % for WSize 40
    [0, 0, 0, 498]   % for WSize 50
];

% Mean and std deviation for SWC
mean1_SWC = {
    [NaN, 41.37, 41.36, 41.36];
    [NaN, NaN, 37.88, 37.87];
    [NaN, NaN, NaN, 34.16];
};
std1_SWC = {
    [NaN, 0.54, 0.6, 0.6];
    [NaN, NaN, 83.58, 0.54];
    [NaN, NaN, NaN, 0.61];
};
mean2_SWC = {
    [NaN, 79.98, 80.39, 80.39];
    [NaN, NaN, 0.54, 83.58];
    [NaN, NaN, NaN, 87.2];
};
std2_SWC = {
    [NaN, 0.13, 0.51, 0.51];
    [NaN, NaN, 0.56, 0.57];
    [NaN, NaN, NaN, 0.45];
};

% Mean and std deviation for HOCo
mean1_HOCo = {
    [NaN, 40, 39.94, 39.94];
    [NaN, NaN, 36.5, 36.55];
    [NaN, NaN, NaN, 33.34 ]
};
std1_HOCo = {
    [NaN, 0, 0.59, 0.59];
    [NaN, NaN, 0.59, 0.56];
    [NaN, NaN, NaN, 0.66]
};
mean2_HOCo = {
    [NaN,80, 81.66, 81.66];
    [NaN, NaN, 84.9, 85.43];
    [NaN, NaN, NaN, 88.73]
};
std2_HOCo = {
    [NaN, 0, 0.59, 0.59];
    [NaN, NaN, 0.59, 0.65];
    [NaN, NaN, NaN, 0.72]
};

% Detection percentage annotations
correct_flip_percentage_SWC = [500, 500, 500]/500*100; % in percentage for three wsizes
correct_flip_percentage_HOCo = [500, 500, 500]/500*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data for True Positives (TP) percentages for SWC and HOCo for transition imposed at 45 and 120 TR
tp_count_SWC = [
    [480, 500, 500, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [0, 0, 500, 500];   % for WSize 40
    [0, 0, 0, 500]      % for WSize 50
];

tp_count_HOCo = [
    [493, 500, 500, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [146, 500, 500, 500];   % for WSize 40
    [0, 500, 500, 500]   % for WSize 50
];

% Mean and std deviation for SWC
mean1_SWC = {
    [46.84, 46.88, 46.88, 46.88];
    [NaN, NaN, 52.12, 52.12 ];
    [NaN, NaN, NaN, 57.41];
};
std1_SWC = {
    [0.37, 0.43, 0.43, 0.43];
    [NaN, NaN, 0.54, 0.54];
    [NaN, NaN, NaN, 0.52];
};
mean2_SWC = {
    [119.09, 119.06, 119.06, 119.06];
    [NaN, NaN, 113.91, 113.91];
    [NaN, NaN, NaN,108.65];
};
std2_SWC = {
    [0.4, 0.43, 0.43, 0.43];
    [NaN, NaN, 0.51, 0.51];
    [NaN, NaN, NaN, 0.56];
};

% Mean and std deviation for HOCo
mean1_HOCo = {
    [46.61, 46.63, 46.63, 46.63];
    [46.97, 47.61, 47.61, 47.61];
    [NaN, 48.97, 48.97, 48.97]
};
std1_HOCo = {
    [0.49, 0.51, 0.51, 0.51];
    [0.16, 0.58, 0.58, 0.58 ];
    [NaN, 0.45, 0.45, 0.45]
};
mean2_HOCo = {
    [118.95, 118.95, 118.95, 118.95];
    [118.03, 117.53, 117.53, 117.53 ];
    [NaN, 115.91, 115.91, 115.91 ]
};
std2_HOCo = {
    [0.36, 0.37, 0.37, 0.37];
    [0.16, 0.55, 0.55, 0.55];
    [NaN, 0.42, 0.42, 0.42]
};

% Detection percentage annotations
correct_flip_percentage_SWC = [500, 500, 500]/500*100; % in percentage for three wsizes
correct_flip_percentage_HOCo = [500, 500, 500]/500*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data for True Positives (TP) percentages for SWC and HOCo for transition imposed at 75 and 120 TR
tp_count_SWC = [
    [406, 500, 500, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [0, 159, 500, 500 ];   % for WSize 40
    [0, 0, 0, 0]      % for WSize 50
];

tp_count_HOCo = [
    [495, 500, 500, 500]; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [14, 500, 500, 500];   % for WSize 40
    [0, 116, 500, 500]   % for WSize 50
];

% Mean and std deviation for SWC
mean1_SWC = {
    [74.84, 74.74, 74.74, 74.74];
    [NaN, 71.12, 70.91, 70.91];
    [NaN, NaN, NaN, NaN];
};
std1_SWC = {
    [0.55, 0.6, 0.6, 0.6];
    [NaN, 0.61, 0.57, 0.57];
    [NaN, NaN, NaN, NaN];
};
mean2_SWC = {
    [121.89, 122.1, 122.1, 122.1];
    [NaN, 124.98, 125.8, 125.8];
    [67.44, 67.44, 67.44, 67.44]; % SWC only found one transition with this mean over 500
};
std2_SWC = {
    [0.32, 0.52, 0.52, 0.52];
    [NaN, 0.14, 0.65, 0.65];
    [0.63, 0.63, 0.63, 0.63]; 
};

% Mean and std deviation for HOCo
mean1_HOCo = {
    [74.58, 74.58, 74.58, 74.58];
    [73, 72.34, 72.34, 72.34];
    [NaN, 70, 69.15, 69.15]
};
std1_HOCo = {
    [0.51, 0.51, 0.51, 0.51];
    [0, 0.52, 0.52, 0.52];
    [NaN, 0, 0.55, 0.55]
};
mean2_HOCo = {
    [121.57, 121.59, 121.59, 121.59];
    [122, 123.07, 123.07, 123.07];
    [NaN, 124.77, 124.7, 124.7]
};
std2_HOCo = {
    [0.5, 0.52, 0.52, 0.52];
    [0, 0.4, 0.4, 0.4];
    [0, 0.42, 0.55, 0.55]
};

% Detection percentage annotations
correct_flip_percentage_SWC = [500, 500, 0]/500*100; % in percentage for three wsizes
correct_flip_percentage_HOCo = [500, 500, 500]/500*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% tp_count_SWC = [
%     [3, 8, 10, 10] ; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
%     [0, 0, 0, 0] ;   % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 40
%     [0, 0, 0, 0]     % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 50
% ];
% 
% tp_count_HOCo = [
%     [0, 33, 285, 421] ; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
%     [0, 0, 16, 59] ;   % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 40
%     [0, 0, 0, 1]       % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 50
% ];
% Data for True Positives (TP) percentages for SWC and HOCo for transition imposed at 15 and 75 TR
% tp_count_SWC = [
%     [59, 117, 123, 123] ; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
%     [0, 0, 3, 3] ;        % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 40
%     [0, 0, 0, 0]          % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 50
% ];
% 
% tp_count_HOCo = [
%     [140, 399, 480, 492] ; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
%     [96, 301, 411, 431] ;  % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 40
%     [27, 187, 365, 410]    % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 50
% ];
% Data for True Positives (TP) percentages for SWC and HOCo for transition imposed at 15 and 120 TR
% tp_count_SWC = [
%     [13, 153, 400, 420] ; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
%     [0, 0, 138, 384] ;    % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 40
%     [0, 0, 0, 131]        % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 50
% ];
% 
% tp_count_HOCo = [
%     [11, 138, 413, 464] ; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
%     [0, 21, 206, 420] ;   % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 40
%     [0, 1, 13, 150]       % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 50
% ];
% Data for True Positives (TP) percentages for SWC and HOCo for transition imposed at 45 and 75 TR
% tp_count_SWC = [
%     [56, 289, 495, 499] ; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
%     [0, 24, 400, 495] ;   % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 40
%     [0, 0, 118, 444]      % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 50
% ];
% 
% tp_count_HOCo = [
%     [1, 75, 431, 488] ; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
%     [0, 1, 132, 443] ;  % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 40
%     [0, 0, 8, 192]      % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 50
% ];
% Data for True Positives (TP) percentages for SWC and HOCo for transition imposed at 45 and 120 TR
% tp_count_SWC = [
%     [239, 467, 500, 500] ; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
%     [8, 108, 460, 500] ;   % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 40
%     [0, 1, 130, 466]      % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 50
% ];
% 
% tp_count_HOCo = [
%     [216, 463, 500, 500] ; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
%     [89, 336, 490, 500] ;  % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 40
%     [7, 201, 493, 500]      % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 50
% ];
% Data for True Positives (TP) percentages for SWC and HOCo for transition imposed at 75 and 120 TR
tp_count_SWC = [
    [236, 456, 500, 500] ; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [52, 256, 480, 481] ;   % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 40
    [0, 50, 91, 94]        % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 50
];

tp_count_HOCo = [
    [226, 462, 495, 496] ; % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 30
    [87, 353, 491, 499] ;  % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 40
    [6, 116, 419, 488]     % Acceptance range: +/-2, +/-5, +/-10, +/-15 TR for WSize 50
];


%% plotting
% Acceptance window ranges
win_ranges = [2, 5, 10, 15]; % +/-2TR, +/-5TR, +/-10TR, +/-15TR

% Window sizes
wsize = [30, 40, 50];

% Plotting
figure;
num_wsize = length(wsize); % Number of window sizes

for i = 1:num_wsize
    subplot(1, num_wsize, i); % Adjust subplots based on number of window sizes
    
    % Smooth line connection for SWC and HOCo using pchip interpolation
    xq = linspace(min(win_ranges), max(win_ranges), 100);
    yq_swc = interp1(win_ranges, tp_count_SWC(i, :), xq, 'pchip'); % pchip interpolation for smooth connection
    yq_hoco = interp1(win_ranges, tp_count_HOCo(i, :), xq, 'pchip');
    
    % Plot the smooth lines for SWC and HOCo
    p1 = plot(xq, yq_swc, '-', 'LineWidth', 1.5, 'Color', 'blue', 'DisplayName', 'SWC (Transition 1)');
    hold on;
    p2 = plot(xq, yq_hoco, '-', 'LineWidth', 1.5, 'Color', 'red', 'DisplayName', 'HOCo (Transition 1)');
    
    
    % Annotation text for SWC and HOCo correct detection percentages
    swc_annotation = sprintf('SWC Correct Detection: %.1f%%', correct_flip_percentage_SWC(i));
    hoco_annotation = sprintf('HOCo Correct Detection: %.1f%%', correct_flip_percentage_HOCo(i));
    
    % Set x-ticks to only show the acceptance window ranges
    set(gca, 'XTick', win_ranges);
    
    % Only show y-axis label on the first subplot
    if i == 1
        ylabel('True Positives (out of 500)', 'FontSize', 12);
    end
    xlabel('Acceptance Window Range (TR)', 'FontSize', 12);
    
    title(['WSize = ', num2str(wsize(i)), ' TR'], 'FontSize', 14);
    grid off;
    ylim([0, 500]); % Set y-axis limits to 0-500 for clarity
    
    % Add a subplot-specific legend with correct detection percentages for each subplot
    lgd = legend({['SWC (Transition 1) (' swc_annotation ')'], ['HOCo (Transition 1) (' hoco_annotation ')'], ...
                  'SWC (Transition 2)', 'HOCo (Transition 2)'}, 'Location', 'southeast');
    lgd.FontSize = 10;
    
    hold off;
end
sgtitle('Transitions at 15 and 120 TR', 'FontSize', 16);

