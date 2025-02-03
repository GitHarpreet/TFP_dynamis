%% Table 4: Mean growth by consecutive years in top 1% based on z_obs

% Calculate growth rates from z_obs
percentile_threshold = 99;  % Define the top 1% threshold percentile

% Initialize array to store the count of years each firm spends in the top 1%
years_in_top_1_percent = zeros(1, N_firms);

% Loop through each time period (1 to T-1)
for t = 1:(T-1)
    % Compute the 99th percentile threshold for z_obs at time t
    top_1_percent_threshold = prctile(z_obs(t, :), percentile_threshold);

    % Identify firms in the top 1% based on z_obs at time t
    firms_in_top_1_percent = z_obs(t, :) >= top_1_percent_threshold;

    % Increment the count of years in top 1% for these firms
    years_in_top_1_percent = years_in_top_1_percent + firms_in_top_1_percent;
end

% Initialize variables to store growth values for firms with consecutive years in the top 1%
growth_1_year = [];
growth_2_years = [];
growth_3_years = [];
growth_4_or_more_years = [];

% Loop through each firm to track consecutive years in the top 1% and calculate growth
for firm = 1:N_firms
    consecutive_count = 0;  % Counter for consecutive years in top 1%

    for t = 1:(T-1)
        if z_obs(t, firm) >= prctile(z_obs(t, :), percentile_threshold)
            % Firm is in the top 1% this year, increment consecutive count
            consecutive_count = consecutive_count + 1;
        else
            % Firm is not in the top 1% this year, record growth and reset counter
            if consecutive_count == 1
                growth_1_year = [growth_1_year; mean(g(max(t-consecutive_count, 1):t-1, firm))];
            elseif consecutive_count == 2
                growth_2_years = [growth_2_years; mean(g(max(t-consecutive_count, 1):t-1, firm))];
            elseif consecutive_count == 3
                growth_3_years = [growth_3_years; mean(g(max(t-consecutive_count, 1):t-1, firm))];
            elseif consecutive_count >= 4
                growth_4_or_more_years = [growth_4_or_more_years; mean(g(max(t-consecutive_count, 1):t-1, firm))];
            end
            consecutive_count = 0;  % Reset count
        end
    end

    % Check if the firm ended with consecutive years in the top 1%
    if consecutive_count == 1
        growth_1_year = [growth_1_year; mean(g(max(T-1-consecutive_count, 1):T-1, firm))];
    elseif consecutive_count == 2
        growth_2_years = [growth_2_years; mean(g(max(T-1-consecutive_count, 1):T-1, firm))];
    elseif consecutive_count == 3
        growth_3_years = [growth_3_years; mean(g(max(T-1-consecutive_count, 1):T-1, firm))];
    elseif consecutive_count >= 4
        growth_4_or_more_years = [growth_4_or_more_years; mean(g(max(T-1-consecutive_count, 1):T-1, firm))];
    end
end

% Calculate mean growth for each category
mean_growth_1_year = mean(growth_1_year);
mean_growth_2_years = mean(growth_2_years);
mean_growth_3_years = mean(growth_3_years);
mean_growth_4_or_more_years = mean(growth_4_or_more_years);

% Create a table to display the results
consecutive_years = {'1 Year'; '2 Years'; '3 Years'; '4 or more years'};
mean_growth_sim = [mean_growth_1_year; mean_growth_2_years; mean_growth_3_years; mean_growth_4_or_more_years];
results_table = table(consecutive_years, mean_growth_sim, 'VariableNames', {'Consecutive_Years', 'Mean_Sim'});

% Display the table
disp('4. Mean growth by consecutive years in Top 1%:');
disp(results_table);
