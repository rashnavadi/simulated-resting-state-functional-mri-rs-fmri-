% function [forward_inflection, backward_inflection, correlation_forward, correlation_backward] = find_inflection_points(X, transition_point)
    % X: data points (n x d matrix) representing time series of
    % correlations, X=[corr(x,y),corr(x,z)corr(y,z)] size: 150x3
    % transitionPoint: the transition point identified by K-means clustering

    figure; plot(all_noisy_x); hold on;plot(all_noisy_y);plot(all_noisy_z)
    title('x, y, z timeseries')
    % hoco beta  values
    figure; plot(dynamic12.bb); hold on; plot(dynamic13.bb); plot(dynamic23.bb)
    title('correlation coefficients (beta values) from HOCo between x, y and z')
    % swc r values
    figure; plot(FC_21); hold on; plot(FC_31); plot(FC_32); 
    title('correlation coefficients (r values) from SWC between x, y and z')
    
    % for HOCo, beta values
    % X = [dynamic12.bb, dynamic13.bb, dynamic23.bb];
    
    % for SWC, r values
    % X = [FC_21, FC_31, FC_32];
    % apply three-point average to remove the noise (as the derivative of the SWC correlations is so jaggedy)
    % Apply three-point moving average to each column
%     X_smoothed = movmean(X, 3);
%     X = X_smoothed;

   % Step 1: Compute the derivative of the correlation
    delta_X = diff(X);
    figure; plot(delta_X(:,1)); hold on; plot(delta_X(:,2)); plot(delta_X(:,3))
    title('HOCo derivative of correlations')
    title('SWC derivative of correlations')


    % Step 2: Find the forward inflection point
    forward_inflection_1 = transitionPoint;
    xy_inflection = delta_X(:,1);
    for j = transitionPoint:size(xy_inflection)
        if j < size(xy_inflection, 1) && (xy_inflection(j) * xy_inflection(j + 1) < 0 || xy_inflection(j) == 0 ) % Change in polarity or zero derivative
            forward_inflection_1 = j + 1;
            break;
        end
    end

    forward_inflection_2 = transitionPoint;
    xz_inflection = delta_X(:, 2);
    for k = transitionPoint:size(delta_X, 1)
        if k < size(xz_inflection, 1) && (xz_inflection(k) * xz_inflection(k + 1) < 0 || xz_inflection(k) == 0) % Change in polarity or zero derivative
            forward_inflection_2 = k + 1;
            break;
        end
    end

    forward_inflection_3 = transitionPoint;
    yz_inflection = delta_X(:, 3);
    for l = transitionPoint:size(delta_X, 1)
        if l < size(yz_inflection, 1) && (yz_inflection(l) * yz_inflection(l + 1) < 0 || yz_inflection(l) == 0 ) % Change in polarity or zero derivative
            forward_inflection_3 = l + 1;
            break;
        end
    end


    % Step 3: Find the backward inflection point
    backward_inflection_1 = transitionPoint;
    for j = transitionPoint-1:-1:2
        if (xy_inflection(j - 1) * xy_inflection(j) < 0 || xy_inflection(j - 1) == 0 ) % Change in polarity or zero derivative
            backward_inflection_1 = j;
            break;
        end
    end

    backward_inflection_2 = transitionPoint;
    for k = transitionPoint-1:-1:2
        if (xz_inflection(k - 1) * xz_inflection(k) < 0 || xz_inflection(k - 1) == 0 ) % Change in polarity or zero derivative
            backward_inflection_2 = k;
            break;
        end
    end

    backward_inflection_3 = transitionPoint;
    for l = transitionPoint-1:-1:2
        if (yz_inflection(l - 1) * yz_inflection(l) < 0 || yz_inflection(l - 1) == 0 ) % Change in polarity or zero derivative
            backward_inflection_3 = l;
            break;
        end
    end


    % Step 4: Calculate the correlation values at the inflection points
    correlation_forward_1 = X(forward_inflection_1, 1);
    correlation_backward_1 = X(backward_inflection_1, 1);

    correlation_forward_2 = X(forward_inflection_2, 2);
    correlation_backward_2 = X(backward_inflection_2, 2);

    correlation_forward_3 = X(forward_inflection_3, 3);
    correlation_backward_3 = X(backward_inflection_3, 3);


    % Display the results
    fprintf('Forward inflection point 1: %d\n', forward_inflection_1);
    fprintf('Backward inflection point 1: %d\n', backward_inflection_1);
    fprintf('Correlation at forward inflection 1: %s\n', mat2str(correlation_forward_1));
    fprintf('Correlation at backward inflection 1: %s\n', mat2str(correlation_backward_1));

    fprintf('Forward inflection point 2: %d\n', forward_inflection_2);
    fprintf('Backward inflection point 2: %d\n', backward_inflection_2);
    fprintf('Correlation at forward inflection 2: %s\n', mat2str(correlation_forward_2));
    fprintf('Correlation at backward inflection 2: %s\n', mat2str(correlation_backward_2));

    fprintf('Forward inflection point 3: %d\n', forward_inflection_3);
    fprintf('Backward inflection point 3: %d\n', backward_inflection_3);
    fprintf('Correlation at forward inflection 3: %s\n', mat2str(correlation_forward_3));
    fprintf('Correlation at backward inflection 3: %s\n', mat2str(correlation_backward_3));



% end
