% Data for Scenario 1 (Single Transition)
% Data for Scenario 1 (Single Transition)
clear 
close all
clc

%% Data for Scenario 1 (Single Transition)
wsize = [30, 40, 50]; % Window sizes
transitions = [15, 45, 75, 120]; % Transition points
swc_scenario1 = [
    24, 29, 34; % Transition at 15 TR
    46, 50, 55; % Transition at 45 TR
    76, 76, 76; % Transition at 75 TR
    117, 111, 106; % Transition at 120 TR
]; % SWC detected transition times
hoco_scenario1 = [
    29, 36, 42; % Transition at 15 TR
    47, 48, 49; % Transition at 45 TR
    76, 76, 76; % Transition at 75 TR
    118, 115, 113; % Transition at 120 TR
];

% Plot Scenario 1
figure;
for i = 1:length(transitions)
    subplot(2, 2, i);
    
    % Plot SWC data points
    p1 = plot(wsize, swc_scenario1(i, :), '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
    hold on;
    
    % Plot HOCo data points
    p2 = plot(wsize, hoco_scenario1(i, :), '-s', 'LineWidth', 1.5, 'MarkerSize', 8);

    % Annotate SWC and HOCo data points with larger font size
   for j = 1:length(wsize)
    if swc_scenario1(i, j) == hoco_scenario1(i, j)
        % If SWC and HOCo are the same, annotate only once above the point
        text(wsize(j), swc_scenario1(i, j) + 3, num2str(swc_scenario1(i, j)), ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 14);
    else
        % SWC annotation (above the point with a small left shift in x direction)
        text(wsize(j) - 0.75, swc_scenario1(i, j) - 5, num2str(swc_scenario1(i, j)), ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 14);
        
        % HOCo annotation (below the point with a small right shift in x direction)
        text(wsize(j) + 0.75, hoco_scenario1(i, j) + 5, num2str(hoco_scenario1(i, j)), ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', 14);
    end
end

    
    hold off;
    
    % Set the X-axis to show only the window sizes 30, 40, 50
    set(gca, 'XTick', wsize);
    xlabel('Window Size (TR)', 'FontSize', 12);
    ylabel('Detected Transition Time (TR)', 'FontSize', 12);
    title(['Transition at ', num2str(transitions(i)), ' TR'], 'FontSize', 14);
    
    % Set the X-ticks to be spaced out for 30, 40, and 50 TR
    xlim([min(wsize)-5, max(wsize)+5]);  % Add some space around the window sizes
    ylim([0, max(max([swc_scenario1(:); hoco_scenario1(:)])) + 10]); % Adjust ylim for visibility
    grid on;
end

% Add a single legend for the entire figure
lgd = legend([p1, p2], 'SWC', 'HOCo', 'FontSize', 12);
lgd.Orientation = 'horizontal';
lgd.Position = [0.3, 0.05, 0.4, 0.05]; % Adjust the position to be below all subplots

% Add a super title for the figure
sgtitle('Scenario 1: Two States and One% Data for Scenario 2 (Two Transitions')


%% two transitions and two states
% Data for Scenario 2 (Two Transitions)
swc_scenario2 = {
    [49, NaN], [51, NaN], [55, NaN];             % Transitions at 15, 45 TR, but SWC only found one transition 
    [76, NaN], [77, NaN], [80, NaN];             % Transitions at 15, 75 TR but SWC only found one transition 
    [22, 114], [26, 109], [31, 104];             % Transitions at 15, 120 TR, SWC found two transitions
    [41, 80], [38, 84], [34, 87];                % Transitions at 45, 75 TR, SWC found two transitions
    [47, 119], [52, 114], [58, 108];             % Transitions at 45, 120 TR
    [75, 122], [71, 126], [67, NaN]              % Transitions at 75, 120 TR, but SWC found only one transition 
};

hoco_scenario2 = {
    [7, 52], [57, NaN], [60, NaN];               % Transitions at 15, 45 TR
    [14, 77], [15, 76], [13, 79];                % Transitions at 15, 75 TR
    [22, 115], [25, 111], [30, 107];             % Transitions at 15, 120 TR
    [40, 82], [36, 86], [34, 89];                % Transitions at 45, 75 TR
    [47, 119], [47, 118], [49, 116];             % Transitions at 45, 120 TR
    [75, 122], [72, 123], [69, 125]              % Transitions at 75, 120 TR
};

% Transition points for the imposed state changes
transitions = {
    [15, 45];
    [15, 75];
    [15, 120];
    [45, 75];
    [45, 120];
    [75, 120];
};

% Window sizes
wsize = [30, 40, 50]; 

% Shades of black for imposed transition lines
transition_colors = [0.1, 0.1, 0.1; 0.6, 0.6, 0.6];

% Plot Scenario 2
figure;
hold on;

% Plot dummy markers for legend
plot(nan, nan, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'DisplayName', 'SWC'); % Dummy for SWC
plot(nan, nan, 's', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'DisplayName', 'HOCo'); % Dummy for HOCo


for i = 1:6 % Six transition combinations in Scenario 2
    subplot(3, 2, i);

    % Draw horizontal lines for the imposed transitions with shades of black
    yline(transitions{i}(1), '--', 'LineWidth', 1.5, 'Color', transition_colors(1,:), 'DisplayName', 'Imposed Transition 1');
    yline(transitions{i}(2), '--', 'LineWidth', 1.5, 'Color', transition_colors(2,:), 'DisplayName', 'Imposed Transition 2');

    hold on;

    for j = 1:length(wsize)
        % SWC may find either one or two transitions, so loop through the
        % two possible detected transitions
        if ~isnan(swc_scenario2{i, j}(1)) % Check if SWC found the first transition
            swc_plot = plot(wsize(j) - 0.5, swc_scenario2{i, j}(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
        end
        if ~isnan(swc_scenario2{i, j}(2)) % Check if SWC found the second transition
            plot(wsize(j) - 0.5, swc_scenario2{i, j}(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
        end
    end

    % Plot HOCo detected transitions (green squares)
    for j = 1:length(wsize)
        if ~isnan(hoco_scenario2{i, j}(1)) % Check if HOCo found the first transition
            hoco_plot = plot(wsize(j) + 0.5, hoco_scenario2{i, j}(1), 's', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
        end
        if ~isnan(hoco_scenario2{i, j}(2)) % Check if HOCo found the second transition
            plot(wsize(j) + 0.5, hoco_scenario2{i, j}(2), 's', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
        end
    end

    % Set the X-axis to show only the window sizes 30, 40, 50
    set(gca, 'XTick', wsize);
    xlabel('Window Size (TR)', 'FontSize', 12);
    ylabel('Detected Transition Time (TR)', 'FontSize', 12);
    title(['Scenario 2: Transition ' num2str(transitions{i}(1)) ' and ' num2str(transitions{i}(2)) ' TR'], 'FontSize', 14);

    % Adjust limits for clarity
    xlim([min(wsize)-5, max(wsize)+5]);
    ylim([0, max([max(cell2mat(swc_scenario2(:))), max(cell2mat(hoco_scenario2(:)))]) + 20]);
    grid on;

    hold off;
end

% Add a single legend for the entire figure
legend([swc_plot, hoco_plot], {'SWC', 'HOCo'}, 'FontSize', 12, 'Location', 'bestoutside');

% Add a super title for the figure
sgtitle('Scenario 2: Two States and Two Transitions (SWC vs. HOCo)', 'FontSize', 16);


