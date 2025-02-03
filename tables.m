z_obs = z_obs_dc;

%% Table 2: Ration of years spent in top 1%
T_years = T - 1;
max_lag = 12;
max_lead = 12;
% store the count of years each firm spends in the top 1%
years_in_top_1_percent = zeros(1, N_firms);

% Loop through each year and calculate the top 1% threshold
for t = 1:T_years
    % Determine the threshold
    top_1_percent_threshold = prctile(z_obs(t, :), 99);

    % firms in the top 1% for this year
    firms_in_top_1_percent = z_obs(t, :) >= top_1_percent_threshold;

    years_in_top_1_percent = years_in_top_1_percent + firms_in_top_1_percent;
end

% ratio of years spent in the top 1% for each firm
ratio_years_in_top_1_percent = years_in_top_1_percent / T_years;

% Count the number of firms that fall within each possible range (0 to 13 years)
count_years_in_top_1_percent = histcounts(years_in_top_1_percent, -0.5:1:13.5);

% Calculate the ratio of firms for each number of years spent in the top 1%
ratio_of_firms = count_years_in_top_1_percent / N_firms;

% Create a table to display the results
yrs = (0:13)'; 
table_2_results = table(yrs, ratio_of_firms', 'VariableNames', {'yrs', 'Number_of_years_Sim'});

% Display the table
disp('2. Ratio of years spent in the top 1% (V):');
disp(table_2_results);

%% Table 3 : Average percentile in the top 1% in the past and the future

average_percentile_past = zeros(max_lag, 1);
average_percentile_future = zeros(max_lead, 1);

% Calculate the thresholds for the top 1% for all years Reference thresholds for lags and leads
top_1_percent_thresholds = prctile(z_obs, 99, 2);
reference_threshold_past = top_1_percent_thresholds(T_years);  
reference_threshold_future = top_1_percent_thresholds(1);     

for lag = 1:max_lag
    total_percentile = 0;  
    count = 0;
    
    % Identify the lagged year
    lag_year = T_years - lag;  % Example: lag = 12 -> year 1; lag = 1 -> year 12
    
    if lag_year > 0
        reference_top_1_firms = find(z_obs(T_years, :) >= reference_threshold_past);
        firm_percentiles = tiedrank(z_obs(lag_year, :)) ./ N_firms;
        total_percentile = sum(firm_percentiles(reference_top_1_firms));
        count = numel(reference_top_1_firms);
    end
    
    % Calculate the average percentile for this lag
    average_percentile_past(lag) = total_percentile / max(count, 1);  % Avoid division by zero
end

for lead = 1:max_lead
    total_percentile = 0;  
    count = 0;             
    
    % Identify the lead year, lead = 12 -> year 13; lead = 1 -> year 2
    lead_year = 1 + lead; 
    
    if lead_year <= T_years
        reference_top_1_firms = find(z_obs(1, :) >= reference_threshold_future);
        firm_percentiles = tiedrank(z_obs(lead_year, :)) ./ N_firms;
        total_percentile = sum(firm_percentiles(reference_top_1_firms));
        count = numel(reference_top_1_firms);
    end
    
    % Calculate the average percentile for this lead and to avoid dividing by 0
    average_percentile_future(lead) = total_percentile / max(count, 1); 
end

% Define lags and leads directly
lags = -max_lag:-1;  
leads = 1:max_lead; 

average_percentile_past = flip(average_percentile_past);

% Combine results
table_lags = table(lags', average_percentile_past, 'VariableNames', {'Years_From_Present', 'Percentile_of_Top_1_Percent_Firms'});

table_leads = table(leads', average_percentile_future, 'VariableNames', {'Years_From_Present', 'Percentile_of_Top_1_Percent_Firms'});

% Combine both tables
table_current = table(0, 1, 'VariableNames', {'Years_From_Present', 'Percentile_of_Top_1_Percent_Firms'}); % Percentile is 1 at reference year
table_3_results = [table_lags; table_current; table_leads];

disp('3. Average percentile of the top 1% in the past and future:');
disp(table_3_results);

%% Table 4: Mean growth by consecutive years in top 1% based on z_obs

% Calculate growth rates from z_obs
g = diff(z_obs); 
percentile_threshold = 99; 

% store the count of years each firm spends in the top 1%
years_in_top_1_percent = zeros(1, N_firms);

% Loop through each time period (1 to T-1)
for t = 1:(T-1)
    % 99th percentile threshold for z_obs at time t
    top_1_percent_threshold = prctile(z_obs(t, :), percentile_threshold);

    % Identify firms in the top 1% based on z_obs
    firms_in_top_1_percent = z_obs(t, :) >= top_1_percent_threshold;

    years_in_top_1_percent = years_in_top_1_percent + firms_in_top_1_percent;
end

% store growth values for firms with consecutive years in the top 1%
growth_1_year = [];
growth_2_years = [];
growth_3_years = [];
growth_4_or_more_years = [];

% Loop through each firm to track consecutive years in the top 1% and calculate growth
for firm = 1:N_firms
    consecutive_count = 0; 

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
table_4_results = table(consecutive_years, mean_growth_sim, 'VariableNames', {'Consecutive_Years', 'Mean_Sim'});

% Display the table
disp('4. Mean growth by consecutive years in Top 1%:');
disp(table_4_results);
%% Exporting the tables
writetable(table_2_results, 'table_2.csv')
writetable(table_3_results, 'table_3.csv')
writetable(table_4_results, 'table_4.csv')