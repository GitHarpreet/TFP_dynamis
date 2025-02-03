clear all;
clc;
close all;

%% Parameters
sigma = [0.10811, 0.27999, 0.09357];  % Standard deviations for each regime
df = 4.0495;  % Degrees of freedom for t-distribution
alpha_z = [0, 0.017038, 0.011278];  % Means for regime 2 and 3, regime 1 is 0
beta_z = [0.7737, 0.64244, 0.95214];  % Persistence for each regime
gamma_z = 0.29729;  % Transition probability adjustment
sigma_u = 0.017714;  % Standard deviation of the measurement error


p1 = 0.0044488;  % Probability of leaving regime 1
p2 = 0.051609;   
p3 = 0.026228;   % Probability of leaving regime 3

% Simulation parameters
N = 11;  % Number of grid points for z_t in each regime
num_regimes = 3;  % Number of regimes
T = 14;  % Firm lifespan in years
burn_in_period = 186;  % Burn-in period
N_firms = 100000;  % Number of firms to simulate


%% Discretize z_t^* using Rouwenhorst method for each regime
z_grids = cell(num_regimes, 1);
P_within = cell(num_regimes, 1);

% Compute the long-run mean for each regime: mu_i = alpha_i / (1 - beta_i)
mu_z = alpha_z ./ (1 - beta_z); 

for i = 1:num_regimes
    [z_grids{i}, P_within{i}] = rouwenhorst(beta_z(i), sigma(i), N, mu_z(i), df);
end

% Transition matrix between regimes (P)
P_between = [1-p1, p1, 0;
             gamma_z*p2, 1-p2, (1-gamma_z)*p2;
             0, p3, 1-p3];

P_between_cumsum = cumsum(P_between, 2); 

% Simulate TFP process with transitions (including 186 burn-in period)
z_obs_star = zeros(T + burn_in_period, N_firms); 
z_obs = zeros(T + burn_in_period, N_firms);  % z_t with measurement error

% Initialize each firm in regime 1
S = ones(1, N_firms);  % Regime state for each firm
for firm = 1:N_firms
    % Random start based on the stationary distribution of the assigned regime
    z_obs_star(1, firm) = mu_z(S(firm)) + trnd(df) * sigma(S(firm));  % Shock drawn from t-distribution
    z_obs(1, firm) = z_obs_star(1, firm) + normrnd(0, sigma_u);  
end

% Simulate TFP process with transitions
for t = 2:(T + burn_in_period)
 
    rand_vals = rand(1, N_firms);  % One random value for each firm

    for firm = 1:N_firms
        prev_S = S(firm);  % Previous regime
        prev_z = z_obs_star(t-1, firm);  % Previous latent TFP value
        
        % Determine the next regime
        next_regime = sum(rand_vals(firm) > P_between_cumsum(prev_S, :)) + 1;
        S(firm) = next_regime;  % Assign the next regime
        epsilon = trnd(df) * sigma(S(firm));  % Regime-specific shock t-dist

  % If the regime changes, use linear interpolation to assign a new grid
       % if prev_S ~= next_regime
            % Interpolate the previous TFP value into the new regime's grid
         %   new_z_grid = z_grids{next_regime};
            
            % Linear interpolation to find the nearest grid point in the new regime
          %  interpolated_z = interp1(z_grids{prev_S}, z_grids{next_regime}, prev_z, 'linear', 'extrap');
            
            % Find the closest grid point in the new regime
          %  [~, nearest_grid_idx] = min(abs(new_z_grid - interpolated_z));
            
            % Assign the TFP value to the nearest grid point in the new regime
           % z_obs_star(t, firm) = new_z_grid(nearest_grid_idx);
       % else
            % Within-regime transition
        %    transition_prob = P_within{S(firm)};
          %  z_obs_star(t, firm) = z_grids{S(firm)}(find(rand < cumsum(transition_prob(round(N/2), :)), 1));
        %end


        % Compute z_star based on the current regime
        if S(firm) == 1
            z_obs_star(t, firm) = beta_z(1) * prev_z + epsilon;
        elseif S(firm) == 2
            z_obs_star(t, firm) = alpha_z(2) + beta_z(2) * prev_z + epsilon;
        elseif S(firm) == 3
            z_obs_star(t, firm) = alpha_z(3) + beta_z(3) * prev_z + epsilon;
        end

        % Add normal measurement error to the observed z_t
        z_obs(t, firm) = z_obs_star(t, firm) + normrnd(0, sigma_u);
    end
end

% Remove the burn-in period

z_obs_star = z_obs_star(burn_in_period+1:end, :);
z_obs = z_obs(burn_in_period+1:end, :);

% Calculate growth as g_t = z_t - z_{t-1} 
g = diff(z_obs);

% Compute the overall mean
overall_mean_growth = mean(g(:));
overall_std_growth = std(g(:));

% Display the results
disp('Overall Mean Growth Rate:');
disp(overall_mean_growth);

disp('Overall Standard Deviation of Growth Rate:');
disp(overall_std_growth);

% percentile
g = diff(z_obs);  % Growth rate for each period and firm
percentile_ranges = [0, 10, 50, 90, 95, 99, 100];

% Initialize storage for yearly mean and variance growth in each percentile range
yearly_mean_growth = zeros(T-1, length(percentile_ranges)-1);
yearly_var_growth = zeros(T-1, length(percentile_ranges)-1);
yearly_skewness_growth = zeros(T-1, length(percentile_ranges)-1);
yearly_kelly_skewness_growth = zeros(T-1, length(percentile_ranges)-1);
yearly_kurtosis_growth = zeros(T-1, length(percentile_ranges)-1);

% Initialize storage for growth rates over the entire period
growth_values_by_percentile = cell(length(percentile_ranges)-1, 1);

% Loop over each time period (we have T-1 periods of growth)
for t = 1:(T-1)
    % Compute the percentile thresholds for z_obs at time t
    percentile_values = prctile(z_obs(t, :), percentile_ranges);

    % For each percentile range, compute mean and variance of growth
    for i = 1:(length(percentile_ranges)-1)
        % Identify firms within the current percentile range
        idx_in_percentile = (z_obs(t, :) >= percentile_values(i)) & (z_obs(t, :) < percentile_values(i+1));
        
        % Extract growth rates for firms in the current percentile range
        growth_in_percentile = g(t, idx_in_percentile);
        
        % Calculate yearly mean and variance for the current percentile range
        if sum(idx_in_percentile) > 0  % Check if there are firms in the bin
            yearly_mean_growth(t, i) = mean(growth_in_percentile);
            yearly_var_growth(t, i) = var(growth_in_percentile);
            yearly_skewness_growth(t, i) = skewness(growth_in_percentile);
            percentiles = prctile(growth_in_percentile, [10, 50, 90]);
            yearly_kelly_skewness_growth(t, i) = ((percentiles(3) - percentiles(2)) - (percentiles(2) - percentiles(1))) / (percentiles(3) - percentiles(1));
            yearly_kurtosis_growth(t, i) = kurtosis(growth_in_percentile);

            % Append growth values for the entire period
            growth_values_by_percentile{i} = [growth_values_by_percentile{i}, growth_in_percentile];
        else
            yearly_mean_growth(t, i) = NaN;  % Handle empty bins
            yearly_var_growth(t, i) = NaN;   % Handle empty bins
            yearly_skewness_growth(t, i) = NaN;
            yearly_kelly_skewness_growth(t, i) = NaN;
            yearly_kurtosis_growth(t, i) = NaN;
        end
    end
end

% Compute the overall mean and variance of growth rates for each percentile range
overall_mean_growth = zeros(length(percentile_ranges)-1, 1);
overall_var_growth = zeros(length(percentile_ranges)-1, 1);
overall_skewness_growth = zeros(length(percentile_ranges)-1, 1);
overall_kelly_skewness_growth = zeros(length(percentile_ranges)-1, 1);
overall_kurtosis_growth = zeros(length(percentile_ranges)-1, 1);

for i = 1:(length(percentile_ranges)-1)
    if ~isempty(growth_values_by_percentile{i})
        overall_mean_growth(i) = mean(growth_values_by_percentile{i});
        overall_var_growth(i) = var(growth_values_by_percentile{i});
        overall_skewness_growth(i) = skewness(growth_values_by_percentile{i});
        percentiles = prctile(growth_values_by_percentile{i}, [10, 50, 90]);
        overall_kelly_skewness_growth(i) = ((percentiles(3) - percentiles(2)) - (percentiles(2) - percentiles(1))) / (percentiles(3) - percentiles(1));
        overall_kurtosis_growth(i) = kurtosis(growth_values_by_percentile{i});
    else
        overall_mean_growth(i) = NaN;
        overall_var_growth(i) = NaN;
        overall_skewness_growth(i) = NaN;
        overall_kelly_skewness_growth(i) = NaN;
        overall_kurtosis_growth(i) = NaN;
    end
end

% Display the yearly results
for t = 1:(T-1)
    fprintf('Year %d:\n', t);
    for i = 1:(length(percentile_ranges)-1)
        fprintf('Percentile range %d%%-%d%%: Mean growth = %f, Variance = %f\n', ...
                percentile_ranges(i), percentile_ranges(i+1), yearly_mean_growth(t, i), yearly_var_growth(t, i));
    end
    disp('--------------------------------');
end

% Display overall results in a table
overall_results_table = table(percentile_ranges(1:end-1)', percentile_ranges(2:end)', ...
                              overall_mean_growth, overall_var_growth, ...
                              overall_skewness_growth, overall_kelly_skewness_growth, ...
                              overall_kurtosis_growth, ...
                              'VariableNames', {'Percentile Start', 'Percentile End', ...
                                                'Mean Growth', 'Variance Growth', ...
                                                'Skewness', 'Kelly Skewness', 'Kurtosis'});

disp('Overall Mean, Variance, Skewness, Kelly Skewness, and Kurtosis of Growth Rates for Each Percentile Range:');
disp(overall_results_table);

writetable(overall_results_table, 'table_1.csv')

% Define the percentile ranges
percentile_labels = {'0-10%', '10-50%', '50-90%', '90-95%', '95-99%', '99-100%'};

mean_tfp_by_percentile = overall_mean_growth;
sd_tfp_by_percentile = sqrt(overall_var_growth);

x_points = 1:length(percentile_labels);  % Equal spacing for each percentile bin

% Plotting a smooth line plot for mean TFP by percentile range
figure;
plot(x_points, mean_tfp_by_percentile, '-o', 'LineWidth', 2, 'MarkerSize', 6);

set(gca, 'XTick', x_points);
xticklabels(percentile_labels); 

% Label and title the plot
xlabel('TFP Percentile');
ylabel('Mean Growth');
title('Mean Growth by TFP Percentile');
grid on;

% Plotting a smooth line plot for mean TFP by percentile range
figure;
plot(x_points, sd_tfp_by_percentile, '-o', 'LineWidth', 2, 'MarkerSize', 6);

set(gca, 'XTick', x_points);
xticklabels(percentile_labels);

% Label and title the plot
xlabel('TFP Percentile');
ylabel('SD Growth');
title('SD by TFP Percentile');
grid on;


%% Rouwenhorst method for discretizing AR(1) process in each regime
function [zgrid, P] = rouwenhorst(rho, sigma_eps, n, mu_z, df)
    % rho: AR(1) persistence parameter
    % sigma_eps: standard deviation of the shocks
    % n: number of grid points for z_t
    % mu_z: long-run mean of z_t (mu_i = alpha_i / (1 - rho))
    % df: degrees of freedom for the t-distribution

    % Calculate q (persistence probability based on rho)
    q = (rho + 1) / 2;  % Probability for persistence

    % Correct sigma_eps for the t-distribution based on df
    correction_factor = sqrt(df / (df - 2));
    sigma_eps_corrected = sigma_eps * correction_factor;

    % Calculate nu based on the unconditional standard deviation of the AR(1) process
    nu = ((n-1) / (1 - rho^2))^(1/2) * sigma_eps_corrected;  % This is the range for z_t

    % Initial transition matrix for two points
    P = [q 1-q; 1-q q];

    % Iteratively construct transition matrix for n points
    for i = 2:n-1
        P = q * [P zeros(i, 1); zeros(1, i+1)] + ...
            (1-q) * [zeros(i, 1) P; zeros(1, i+1)] + ...
            (1-q) * [zeros(1, i+1); P zeros(i, 1)] + ...
            q * [zeros(1, i+1); zeros(i, 1) P];
        P(2:i,:) = P(2:i,:) / 2;  % Normalize
    end

    % Create the grid for z centered around the long-run mean mu_z
    zgrid = linspace(mu_z - nu, mu_z + nu, n);  % Use nu to cover the range of possible values
end
