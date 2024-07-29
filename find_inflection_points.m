function [forward_inflection, backward_inflection, correlation_forward, correlation_backward] = find_inflection_points(X, transition_point)
    % X: data points (n x d matrix) representing time series of
    % correlations, X=[x,y,z] size: 150x3
    % transition_point: the transition point identified by K-means clustering

    % Step 1: Compute the derivative of the correlation
    delta_X = diff(X);

    % Step 2: Find the forward inflection point
    forward_inflection = transition_point;
    for j = transition_point:size(delta_X, 1)
        if j < size(delta_X, 1) && (delta_X(j) * delta_X(j + 1) < 0 || delta_X(j) == 0 || delta_X(j + 1) == 0) % Change in polarity or zero derivative
            forward_inflection = j + 1;
            break;
        end
    end

    % Step 3: Find the backward inflection point
    backward_inflection = transition_point;
    for j = transition_point:-1:2
        if (delta_X(j - 1) * delta_X(j) < 0 || delta_X(j - 1) == 0 || delta_X(j) == 0) % Change in polarity or zero derivative
            backward_inflection = j;
            break;
        end
    end

    % Step 4: Calculate the correlation values at the inflection points
    correlation_forward = X(forward_inflection, :);
    correlation_backward = X(backward_inflection, :);

    % Display the results
    fprintf('Forward inflection point: %d\n', forward_inflection);
    fprintf('Backward inflection point: %d\n', backward_inflection);
    fprintf('Correlation at forward inflection: %s\n', mat2str(correlation_forward));
    fprintf('Correlation at backward inflection: %s\n', mat2str(correlation_backward));
end
