% Define percentile ranges
percentile_ranges = [0, 10, 50, 90, 95, 99, 100];

% Number of percentile bins
num_bins = length(percentile_ranges) - 1;

% Initialize the 5-year transition matrix
transition_matrix_5yr = zeros(num_bins, num_bins);

% Calculate 5-year growth rates (difference every 5 years)
z_obs_5yr = z_obs(1:5:end, :);  % Subsample every 5th year
num_periods_5yr = size(z_obs_5yr, 1) - 1;  % Number of 5-year intervals

% Loop over each 5-year interval
for t = 1:num_periods_5yr
    % Compute the percentiles for the current and next period
    percentiles_t = prctile(z_obs_5yr(t, :), percentile_ranges);
    percentiles_tplus5 = prctile(z_obs_5yr(t+1, :), percentile_ranges);

    % Determine the bin for each firm in the current and next period
    bins_t = discretize(z_obs_5yr(t, :), percentiles_t);
    bins_tplus5 = discretize(z_obs_5yr(t+1, :), percentiles_tplus5);

    % Accumulate transitions between bins
    for i = 1:num_bins
        for j = 1:num_bins
            % Count transitions from bin i to bin j
            transition_matrix_5yr(i, j) = transition_matrix_5yr(i, j) + ...
                sum(bins_t == i & bins_tplus5 == j);
        end
    end
end

% Normalize each row to get probabilities
transition_matrix_5yr = transition_matrix_5yr ./ sum(transition_matrix_5yr, 2);

% Display the transition matrix
disp('5-Year Transition Matrix for TFP (Percentile Bins):');
disp(array2table(transition_matrix_5yr, ...
    'VariableNames', strcat('To_', string(percentile_ranges(2:end))), ...
    'RowNames', strcat('From_', string(percentile_ranges(1:end-1)))));





% Define percentile ranges
percentile_ranges = [0, 10, 50, 90, 95, 99, 100];

% Number of percentile bins
num_bins = length(percentile_ranges) - 1;

% Initialize the 1-year transition matrix
transition_matrix_1yr = zeros(num_bins, num_bins);

% Calculate 1-year growth rates (difference every 1 years)
num_periods_1yr = size(z_obs, 1) - 1;  % Number of 1-year intervals

% Loop over each 5-year interval
for t = 1:num_periods_1yr
    % Compute the percentiles for the current and next period
    percentiles_t = prctile(z_obs(t, :), percentile_ranges);
    percentiles_tplus1 = prctile(z_obs(t+1, :), percentile_ranges);

    % Determine the bin for each firm in the current and next period
    bins_t = discretize(z_obs(t, :), percentiles_t);
    bins_tplus1 = discretize(z_obs(t+1, :), percentiles_tplus1);

    % Accumulate transitions between bins
    for i = 1:num_bins
        for j = 1:num_bins
            % Count transitions from bin i to bin j
            transition_matrix_1yr(i, j) = transition_matrix_1yr(i, j) + ...
                sum(bins_t == i & bins_tplus1 == j);
        end
    end
end

% Normalize each row to get probabilities
transition_matrix_1yr = transition_matrix_1yr ./ sum(transition_matrix_1yr, 2);

% Display the transition matrix
disp('1-Year Transition Matrix for TFP (Percentile Bins):');
disp(array2table(transition_matrix_1yr, ...
    'VariableNames', strcat('To_', string(percentile_ranges(2:end))), ...
    'RowNames', strcat('From_', string(percentile_ranges(1:end-1)))));

