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

sigma = sigma * sqrt((df/(df-2)));

p1 = 0.0044488;  % Probability of leaving regime 1
p2 = 0.051609;   
p3 = 0.026228;   % Probability of leaving regime 3  

% Simulation parameters
N = 15;  % Number of grid points for z_t in each regime
num_regimes = 3;  % Number of regimes
T = 14;  % Firm lifespan in years
burn_in_period = 186;  % Burn-in period
N_firms = 100000;  % Number of firms to simulate

z_grids = cell(num_regimes, 1);
P_within = cell(num_regimes, 1);

% Compute the long-run mean for each regime: mu_i = alpha_i / (1 - beta_i)
mu_z = alpha_z ./ (1 - beta_z); 

%for i = 1:num_regimes
%    [z_grids{i}, P_within{i}] = tauchen(beta_z(i), sigma(i), N, mu_z(i), df);
%end
for i  = 1: num_regimes

[P_within{i},z_grids{i}] = discreteAR(mu_z(i),beta_z(i),sigma(i),15,'even',2,2);
end

% Transition matrix between regimes (P)
P_between = [1-p1, p1, 0;
             gamma_z*p2, 1-p2, (1-gamma_z)*p2;
             0, p3, 1-p3];

P_between_cumsum = cumsum(P_between, 2); 

% Simulate TFP process with transitions (including 186 burn-in period)
z_obs_star = zeros(T + burn_in_period, N_firms); 
z_obs_dc = zeros(T + burn_in_period, N_firms);  % z_t with measurement error

% Initialize each firm in regime 1
S = ones(1, N_firms);  % Regime state for each firm

for firm = 1:N_firms
    % Random start based on the stationary distribution of the assigned regime
    z_obs_star(1, firm) = mu_z(S(firm)) + trnd(df) * sigma(S(firm));  % Shock drawn from t-distribution
    z_obs_dc(1, firm) = z_obs_star(1, firm) + normrnd(0, sigma_u);  
end

start_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
disp(['Simulation started at: ', char(start_time)]);

% Start the CPU timer
cpu_start = tic;


% Simulate TFP process with transitions
for t = 2:(T + burn_in_period)
    fprintf('Current time step: %d\n', t);  % Print time step for debugging
    rand_vals = rand(1, N_firms);  % Random values for regime transitions

    for firm = 1:N_firms
        prev_S = S(firm);  % Previous regime
        prev_z = z_obs_star(t-1, firm);  % Previous z_t^*

        % Determine the next regime
        next_regime = sum(rand_vals(firm) > P_between_cumsum(prev_S, :)) + 1;
        S(firm) = next_regime;

        % Get grids for current and next regimes
        current_grid = z_grids{prev_S};
        next_grid = z_grids{next_regime};
        epsilon = trnd(df) * sigma(S(firm));  % Regime-specific shock t-dist

        % Handle regime transition
% Handle regime transition
if prev_S ~= next_regime
    % Find the two closest grid points in the next regime
    [~, idx] = sort(abs(next_grid - prev_z), 'ascend');
    closest_points = idx(1:2);  % Indices of the two closest grid points
    
    % Get the values of the two closest grid points
    point_a = next_grid(closest_points(1));  % First closest point
    point_b = next_grid(closest_points(2));  % Second closest point
    
    % Compute probabilities based on relative distances
    prob_a = (point_b - prev_z) / (point_b - point_a);  % Probability for point_a
    prob_b = (prev_z - point_a) / (point_b - point_a);  % Probability for point_b
    
    % Normalize probabilities (optional, since they already sum to 1)
    prob = [prob_a, prob_b];
    cumulative_probs = cumsum(prob);  % Cumulative probabilities


    % Assign the next grid point based on random value
    random_value = rand;
    if random_value <= cumulative_probs(1)
        selected_idx = closest_points(1);
    else
        selected_idx = closest_points(2);
    end

    % Assign the new z_t^* based on selected index
    z_obs_star(t, firm) = next_grid(selected_idx); % + epsilon;
else
    % Within-regime transition: use existing logic
    [~, current_idx] = min(abs(current_grid - prev_z));  % Closest grid point in current regime
    transition_probs = P_within{prev_S}(current_idx, :);
    cumulative_probs = cumsum(transition_probs);
    random_value = rand;

    % Direct indexing for within-regime transitions
    for i = 1:length(cumulative_probs)
        if random_value <= cumulative_probs(i)
            next_z_idx = i;  % Selected index
            break;
        end
    end

    % Assign the new z_t^* within the same regime
    z_obs_star(t, firm) = current_grid(next_z_idx); % + epsilon;
end


        % Add measurement error
        z_obs_dc(t, firm) = z_obs_star(t, firm) + normrnd(0, sigma_u);
    end
end

% Stop the CPU timer
cpu_elapsed = toc(cpu_start);

% Display the end time
end_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
disp(['Simulation ended at: ', char(end_time)]);

% Display the CPU elapsed time
disp(['CPU time elapsed: ', num2str(cpu_elapsed), ' seconds']);



% Compute the overall mean

g = diff(z_obs_dc);  % Growth rate for each period and firm
percentile_ranges = [0, 10, 50, 90, 95, 99, 100];

% percentile
g = diff(z_obs_dc);  % Growth rate for each period and firm
percentile_ranges = [0, 10, 50, 90, 95, 99, 100];

% Initialize storage for yearly mean and variance growth in each percentile range
yearly_mean_growth = zeros(T-1, length(percentile_ranges)-1);
yearly_sd_growth = zeros(T-1, length(percentile_ranges)-1);
yearly_skewness_growth = zeros(T-1, length(percentile_ranges)-1);
yearly_kelly_skewness_growth = zeros(T-1, length(percentile_ranges)-1);
yearly_kurtosis_growth = zeros(T-1, length(percentile_ranges)-1); 

% Initialize storage for growth rates over the entire period
growth_values_by_percentile = cell(length(percentile_ranges)-1, 1);

% Loop over each time period (we have T-1 periods of growth)
for t = 1:(T-1)
    % Compute the percentile thresholds for z_obs at time t
    percentile_values = prctile(z_obs_dc(t, :), percentile_ranges);

    % For each percentile range, compute mean and variance of growth
    for i = 1:(length(percentile_ranges)-1)
        % Identify firms within the current percentile range
        idx_in_percentile = (z_obs_dc(t, :) >= percentile_values(i)) & (z_obs_dc(t, :) < percentile_values(i+1));
        
        % Extract growth rates for firms in the current percentile range
        growth_in_percentile = g(t, idx_in_percentile);
        
        % Calculate yearly mean and variance for the current percentile range
        if sum(idx_in_percentile) > 0  % Check if there are firms in the bin
            yearly_mean_growth(t, i) = mean(growth_in_percentile);
            yearly_sd_growth(t, i) = std(growth_in_percentile);
            yearly_skewness_growth(t, i) = skewness(growth_in_percentile);
            percentiles = prctile(growth_in_percentile, [10, 50, 90]);
            yearly_kelly_skewness_growth(t, i) = ((percentiles(3) - percentiles(2)) - (percentiles(2) - percentiles(1))) / (percentiles(3) - percentiles(1));
            yearly_kurtosis_growth(t, i) = kurtosis(growth_in_percentile);

            % Append growth values for the entire period
            growth_values_by_percentile{i} = [growth_values_by_percentile{i}, growth_in_percentile];
        else
            yearly_mean_growth(t, i) = NaN;  % Handle empty bins
            yearly_sd_growth(t, i) = NaN;   % Handle empty bins
            yearly_skewness_growth(t, i) = NaN;
            yearly_kelly_skewness_growth(t, i) = NaN;
            yearly_kurtosis_growth(t, i) = NaN;
        end
    end
end

% Compute the overall mean and variance of growth rates for each percentile range
overall_mean_growth = zeros(length(percentile_ranges)-1, 1);
overall_std_growth = zeros(length(percentile_ranges)-1, 1);
overall_skewness_growth = zeros(length(percentile_ranges)-1, 1);
overall_kelly_skewness_growth = zeros(length(percentile_ranges)-1, 1);
overall_kurtosis_growth = zeros(length(percentile_ranges)-1, 1);

for i = 1:(length(percentile_ranges)-1)
    if ~isempty(growth_values_by_percentile{i})
        overall_mean_growth(i) = mean(growth_values_by_percentile{i});
        overall_std_growth(i) = std(growth_values_by_percentile{i});
        overall_skewness_growth(i) = skewness(growth_values_by_percentile{i});
        percentiles = prctile(growth_values_by_percentile{i}, [10, 50, 90]);
        overall_kelly_skewness_growth(i) = ((percentiles(3) - percentiles(2)) - (percentiles(2) - percentiles(1))) / (percentiles(3) - percentiles(1));
        overall_kurtosis_growth(i) = kurtosis(growth_values_by_percentile{i});
    else
        overall_mean_growth(i) = NaN;
        overall_std_growth(i) = NaN;
        overall_skewness_growth(i) = NaN;
        overall_kelly_skewness_growth(i) = NaN;
        overall_kurtosis_growth(i) = NaN;
    end
end

% Display the yearly results
for t = 1:(T-1)
    fprintf('Year %d:\n', t);
    for i = 1:(length(percentile_ranges)-1)
        fprintf('Percentile range %d%%-%d%%: Mean growth = %f, std = %f\n', ...
                percentile_ranges(i), percentile_ranges(i+1), yearly_mean_growth(t, i), yearly_sd_growth(t, i));
    end
    disp('--------------------------------');
end

% Display overall results in a table
overall_results_table = table(percentile_ranges(1:end-1)', percentile_ranges(2:end)', ...
                              overall_mean_growth, overall_std_growth, ...
                              overall_skewness_growth, overall_kelly_skewness_growth, ...
                              overall_kurtosis_growth, ...
                              'VariableNames', {'Percentile Start', 'Percentile End', ...
                                                'Mean Growth', 'Std Growth', ...
                                                'Skewness', 'Kelly Skewness', 'Kurtosis'});

disp('Overall Mean, Std, Skewness, Kelly Skewness, and Kurtosis of Growth Rates for Each Percentile Range:');
disp(overall_results_table);


% Define the percentile ranges
percentile_labels = {'0-10%', '10-50%', '50-90%', '90-95%', '95-99%', '99-100%'};

mean_tfp_by_percentile = overall_mean_growth;
sd_tfp_by_percentile = sqrt(overall_std_growth);

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


% Initialize storage for yearly statistics of the top 0.1%
yearly_mean_top_01 = zeros(T-1, 1);
yearly_sd_top_01 = zeros(T-1, 1);
yearly_skewness_top_01 = zeros(T-1, 1);
yearly_kelly_skewness_top_01 = zeros(T-1, 1);
yearly_kurtosis_top_01 = zeros(T-1, 1);

% Storage for all growth rates in top 0.1% across years
growth_values_top_01 = [];

% Loop over each time period (T-1)
for t = 1:(T-1)
    % Find the threshold for the top 0.1% firms at time t
    top_01_threshold = prctile(z_obs_dc(t, :), 99.9); 
    
    % Identify firms that are in the top 0.1% at time t
    idx_top_01 = (z_obs_dc(t, :) >= top_01_threshold);
    
    % Extract growth rates for these firms
    growth_top_01 = g(t, idx_top_01);
    
    % Compute statistics if there are firms in the top 0.1%
    if sum(idx_top_01) > 0
        yearly_mean_top_01(t) = mean(growth_top_01);
        yearly_sd_top_01(t) = std(growth_top_01);
        yearly_skewness_top_01(t) = skewness(growth_top_01);
        percentiles = prctile(growth_top_01, [10, 50, 90]);
        yearly_kelly_skewness_top_01(t) = ((percentiles(3) - percentiles(2)) - (percentiles(2) - percentiles(1))) / (percentiles(3) - percentiles(1));
        yearly_kurtosis_top_01(t) = kurtosis(growth_top_01);

        % Store growth values for the entire period
        growth_values_top_01 = [growth_values_top_01, growth_top_01];
    else
        yearly_mean_top_01(t) = NaN;
        yearly_sd_top_01(t) = NaN;
        yearly_skewness_top_01(t) = NaN;
        yearly_kelly_skewness_top_01(t) = NaN;
        yearly_kurtosis_top_01(t) = NaN;
    end
end

% Compute overall statistics for the top 0.1% firms over all periods
if ~isempty(growth_values_top_01)
    overall_mean_top_01 = mean(growth_values_top_01);
    overall_sd_top_01 = std(growth_values_top_01);
    overall_skewness_top_01 = skewness(growth_values_top_01);
    percentiles = prctile(growth_values_top_01, [10, 50, 90]);
    overall_kelly_skewness_top_01 = ((percentiles(3) - percentiles(2)) - (percentiles(2) - percentiles(1))) / (percentiles(3) - percentiles(1));
    overall_kurtosis_top_01 = kurtosis(growth_values_top_01);
else
    overall_mean_top_01 = NaN;
    overall_sd_top_01 = NaN;
    overall_skewness_top_01 = NaN;
    overall_kelly_skewness_top_01 = NaN;
    overall_kurtosis_top_01 = NaN;
end

% Display the yearly results
fprintf('\nYearly Growth Statistics for the Top 0.1%% (99.9-100%% percentile)\n');
fprintf('------------------------------------------------------------\n');
fprintf('Year | Mean Growth | Std Dev | Skewness | Kelly Skewness | Kurtosis\n');
fprintf('------------------------------------------------------------\n');

for t = 1:(T-1)
    fprintf('%4d | %10.6f | %8.6f | %8.6f | %10.6f | %8.6f\n', ...
        t, yearly_mean_top_01(t), yearly_sd_top_01(t), ...
        yearly_skewness_top_01(t), yearly_kelly_skewness_top_01(t), ...
        yearly_kurtosis_top_01(t));
end

% Display overall statistics
fprintf('\nOverall Growth Statistics for the Top 0.1%% (99.9-100%% percentile)\n');
fprintf('------------------------------------------------------------\n');
fprintf('Mean Growth    : %10.6f\n', overall_mean_top_01);
fprintf('Std Deviation  : %10.6f\n', overall_sd_top_01);
fprintf('Skewness       : %10.6f\n', overall_skewness_top_01);
fprintf('Kelly Skewness : %10.6f\n', overall_kelly_skewness_top_01);
fprintf('Kurtosis       : %10.6f\n', overall_kurtosis_top_01);


%%
% Initialize top 1% threshold (z_99)
z_99_current = prctile(z_obs_dc(:), 99); % Initial guess
tol = 1e-6;  % Convergence tolerance
max_iter = 100;  % Maximum iterations
z_99_prev = 0;  % Previous threshold
iter = 0;

% Iterative adjustment for top 1% threshold
while abs(z_99_current - z_99_prev) > tol && iter < max_iter
    iter = iter + 1;
    z_99_prev = z_99_current;
    
    % Identify firms in the top 1% based on current threshold
    idx_top_1_percent = z_obs_dc(1:end-1,:) > z_99_current;
    
    % Compute growth rates for top 1%
    growth_top_1_percent = g(idx_top_1_percent);
    
    % Update the threshold based on new data
    z_99_current = prctile(z_obs_dc(idx_top_1_percent), 99);
    
    % Print iteration details
    fprintf('Iteration: %d, z_99: %.6f, Mean Growth: %.6f\n', ...
            iter, z_99_current, mean(growth_top_1_percent, 'omitnan'));
end

if iter == max_iter
    warning('Maximum iterations reached. Convergence may not have been achieved.');
else
    fprintf('Converged after %d iterations. Final z_99: %.6f\n', iter, z_99_current);
end

%%
function [zgrid, P] = tauchen(rho, sigma_eps, n, mu_z, df)
    % rho: AR(1) persistence parameter
    % sigma_eps: standard deviation of the shocks
    % n: number of grid points for z_t
    % mu_z: long-run mean of z_t (mu_i = alpha_i / (1 - rho))
    % df: degrees of freedom for the t-distribution

    % Correct sigma_eps for the t-distribution based on df
    correction_factor = sqrt(df / (df - 2));
    sigma_eps_corrected = sigma_eps * correction_factor;

    % Calculate the standard deviation of the stationary AR(1) process
    sigma_z = sigma_eps_corrected / sqrt(1 - rho^2);

    % Create grid points symmetrically around mu_z
    z_min = mu_z - 3 * sigma_z;  % Lower bound
    z_max = mu_z + 3 * sigma_z;  % Upper bound
    zgrid = linspace(z_min, z_max, n);

    % Step size between grid points
    delta = zgrid(2) - zgrid(1);

    % Initialize transition matrix
    P = zeros(n, n);

    % Compute transition probabilities
    for i = 1:n
        for j = 1:n
            if j == 1
                % Lower bound probability
                P(i, j) = tcdf((zgrid(j) - rho * zgrid(i) + delta / 2) / sigma_eps_corrected, df);
            elseif j == n
                % Upper bound probability
                P(i, j) = 1 - tcdf((zgrid(j) - rho * zgrid(i) - delta / 2) / sigma_eps_corrected, df);
            else
                % Middle probabilities
                P(i, j) = tcdf((zgrid(j) - rho * zgrid(i) + delta / 2) / sigma_eps_corrected, df) - ...
                          tcdf((zgrid(j) - rho * zgrid(i) - delta / 2) / sigma_eps_corrected, df);
            end
        end
    end
end
%%



function [P,X] = discreteAR(mu,rho,sigma,Nm,method,nMoments,nSigmas)

% discretize AR(1) process with Gaussian shocks and various grids
% no need to use this because it is a special case of discreteGMAR.m
% define conditional central moments
T1 = 0;
T2 = sigma^2;
T3 = 0;
T4 = 3*sigma^4;

TBar = [T1 T2 T3 T4]'; % vector of conditional central moments

% Default number of moments to match is 2
if nargin == 4
    nMoments = 2;
end

% define grid spacing parameter if not provided
if nargin < 7
    if abs(rho) <= 1-2/(Nm-1)
        nSigmas = sqrt(2*(Nm-1));
    else
        nSigmas = sqrt(Nm-1);
    end
end

% Check that Nm is a valid number of grid points
if ~isnumeric(Nm) || Nm < 3 || rem(Nm,1) ~= 0
    error('Nm must be a positive integer greater than 3')
end

% Check that nMoments is a valid number
if ~isnumeric(nMoments) || nMoments < 1 || nMoments > 4 || ~((rem(nMoments,1) == 0) || (nMoments == 1))
    error('nMoments must be either 1, 2, 3, 4')
end

sigmaX = sigma/sqrt(1-rho^2); % unconditional standard deviation

switch method
    case 'even'
        X = linspace(mu-nSigmas*sigmaX,mu+nSigmas*sigmaX,Nm);
        W = ones(1,Nm);
    case 'gauss-legendre'
        [X,W] = legpts(Nm,[mu-nSigmas*sigmaX,mu+nSigmas*sigmaX]);
        X = X';
    case 'clenshaw-curtis'
        [X,W] = fclencurt(Nm,mu-nSigmas*sigmaX,mu+nSigmas*sigmaX);
        X = fliplr(X');
        W = fliplr(W');
    case 'gauss-hermite'
        [X,W] = GaussHermite(Nm);
        X = mu+sqrt(2)*sigma*X';
        W = W'./sqrt(pi);
end

P = NaN(Nm);
scalingFactor = max(abs(X));
kappa = 1e-8;

for ii = 1:Nm
    
    condMean = mu*(1-rho)+rho*X(ii); % conditional mean
    switch method % define prior probabilities
        case 'gauss-hermite'
            q = W;
        otherwise
            q = W.*normpdf(X,condMean,sigma);
    end
    
    if any(q < kappa)
        q(q < kappa) = kappa; % replace by small number for numerical stability
    end

if nMoments == 1 % match only 1 moment
            P(ii,:) = discreteApproximation(X,@(x)(x-condMean)/scalingFactor,TBar(1)./scalingFactor,q,0);
else % match 2 moments first
    [p,lambda,momentError] = discreteApproximation(X,@(x) [(x-condMean)./scalingFactor;...
        ((x-condMean)./scalingFactor).^2],...
        TBar(1:2)./(scalingFactor.^(1:2)'),q,zeros(2,1));
    if norm(momentError) > 1e-5 % if 2 moments fail, then just match 1 moment
                warning('Failed to match first 2 moments. Just matching 1.')
                P(ii,:) = discreteApproximation(X,@(x)(x-condMean)/scalingFactor,0,q,0);
    elseif nMoments == 2
        P(ii,:) = p;
    elseif nMoments == 3 % 3 moments
    [pnew,~,momentError] = discreteApproximation(X,@(x) [(x-condMean)./scalingFactor;...
        ((x-condMean)./scalingFactor).^2;((x-condMean)./scalingFactor).^3],...
        TBar(1:3)./(scalingFactor.^(1:3)'),q,[lambda;0]);
    if norm(momentError) > 1e-5
        warning('Failed to match first 3 moments.  Just matching 2.')
        P(ii,:) = p;
    else P(ii,:) = pnew;
    end
    else % 4 moments
    [pnew,~,momentError] = discreteApproximation(X,@(x) [(x-condMean)./scalingFactor;...
        ((x-condMean)./scalingFactor).^2; ((x-condMean)./scalingFactor).^3;...
        ((x-condMean)./scalingFactor).^4],TBar./(scalingFactor.^(1:4)'),q,[lambda;0;0]);
    if norm(momentError) > 1e-5
        %warning('Failed to match first 4 moments.  Just matching 3.')
        [pnew,~,momentError] = discreteApproximation(X,@(x) [(x-condMean)./scalingFactor;...
        ((x-condMean)./scalingFactor).^2;((x-condMean)./scalingFactor).^3],...
        TBar(1:3)./(scalingFactor.^(1:3)'),q,[lambda;0]);
    if norm(momentError) > 1e-5
        warning('Failed to match first 3 moments.  Just matching 2.')
        P(ii,:) = p;
    else P(ii,:) = pnew;
        warning('Failed to match first 4 moments.  Just matching 3.')
    end
    else P(ii,:) = pnew;
    end
    end
end

end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discreteApproximation
% (c) 2016 Leland E. Farmer and Alexis Akira Toda
% 
% Purpose: 
%       Compute a discrete state approximation to a distribution with known
%       moments, using the maximum entropy procedure proposed in Tanaka and
%       Toda (2013)
%
% Usage:
%       [p,lambdaBar,momentError] = discreteApproximation(D,T,TBar,q,lambda0)
%
% Inputs:
% D         - (K x N) matrix of grid points. K is the dimension of the
%             domain. N is the number of points at which an approximation
%             is to be constructed.
% T         - A function handle which should accept arguments of dimension
%             (K x N) and return an (L x N) matrix of moments evaluated at
%             each grid point, where L is the number of moments to be
%             matched.
% TBar      - (L x 1) vector of moments of the underlying distribution
%             which should be matched
% Optional:
% q         - (1 X N) vector of prior weights for each point in D. The
%             default is for each point to have an equal weight.
% lambda0   - (L x 1) vector of initial guesses for the dual problem
%             variables. The default is a vector of zeros.
%
% Outputs:
% p         - (1 x N) vector of probabilties assigned to each grid point in
%             D.
% lambdaBar - (L x 1) vector of dual problem variables which solve the
%             maximum entropy problem
% momentError - (L x 1) vector of errors in moments (defined by moments of
%               discretization minus actual moments)
%
% Version 1.2: June 7, 2016
%
% Version 1.3: May 26, 2019
%
% Changed algorithm to 'trust-region' to use Hessian
%
% Version 1.4: September 27, 2023
%
% Changed fminunc option for Matlab 2023b
% Display warning if moment error is large
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%
function [p,lambdaBar,momentError] = discreteApproximation(D,T,TBar,q,lambda0)

% Input error checking

if nargin < 3
    error('You must provide at least 3 arguments to discreteApproximation.')
end

N = size(D,2);

Tx = T(D);
L = size(Tx,1);

if size(Tx,2) ~= N || length(TBar) ~= L
    error('Dimension mismatch')
end

% Set default parameters if not provided
if nargin < 4
    q = ones(1,N)./N; % uniform distribution
end
if nargin < 5
    lambda0 = zeros(L,1);
end

% Compute maximum entropy discrete distribution

%options = optimset('TolFun',1e-10,'TolX',1e-10,'Display','off','GradObj','on','Hessian','on');
%options = optimset('TolFun',1e-10,'TolX',1e-10,'Display','off','Algorithm','trust-region');
options = optimoptions('fminunc','TolFun',1e-10,'TolX',1e-10,'Display','off',...
    'Algorithm','trust-region','SpecifyObjectiveGradient',true);

% Sometimes the algorithm fails to converge if the initial guess is too far
% away from the truth. If this occurs, the program tries an initial guess
% of all zeros.
try
    lambdaBar = fminunc(@(lambda) entropyObjective(lambda,Tx,TBar,q),lambda0,options);
catch
    warning('Failed to find a solution from provided initial guess. Trying new initial guess.')
    lambdaBar = fminunc(@(lambda) entropyObjective(lambda,Tx,TBar,q),zeros(size(lambda0)),options);
end

% Compute final probability weights and moment errors
[obj,gradObj] = entropyObjective(lambdaBar,Tx,TBar,q);
Tdiff = Tx-repmat(TBar,1,N);
p = (q.*exp(lambdaBar'*Tdiff))./obj;
momentError = gradObj./obj;

if norm(momentError) > 1e-5
    warning('Large moment error. Consider increasing number of points or expanding domain')

end
end

%% entropyObjective
% (c) 2016 Leland E. Farmer and Alexis Akira Toda
% 
% Purpose: 
%       Compute the maximum entropy objective function used in
%       discreteApproximation
%
% Usage:
%       obj = entropyObjective(lambda,Tx,TBar,q)
%
% Inputs:
% lambda    - (L x 1) vector of values of the dual problem variables
% Tx        - (L x N) matrix of moments evaluated at the grid points
%             specified in discreteApproximation
% TBar      - (L x 1) vector of moments of the underlying distribution
%             which should be matched 
% q         - (1 X N) vector of prior weights for each point in the grid.
%
% Outputs:
% obj       - scalar value of objective function evaluated at lambda
% Optional (useful for optimization routines):
% gradObj   - (L x 1) gradient vector of the objective function evaluated
%             at lambda
% hessianObj- (L x L) hessian matrix of the objective function evaluated at
%             lambda
%
% Version 1.2: June 7, 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [obj,gradObj,hessianObj] = entropyObjective(lambda,Tx,TBar,q)

% Some error checking

if nargin < 4
    error('You must provide 4 arguments to entropyObjective.')
end

[L,N] = size(Tx);

if length(lambda) ~= L || length(TBar) ~= L || length(q) ~= N
    error('Dimensions of inputs are not compatible.')
end

% Compute objective function

Tdiff = Tx-repmat(TBar,1,N);
temp = q.*exp(lambda'*Tdiff);
obj = sum(temp);

% Compute gradient of objective function

if nargout > 1
    temp2 = bsxfun(@times,temp,Tdiff);
    gradObj = sum(temp2,2);
end

% Compute hessian of objective function

if nargout > 2
    hessianObj = temp2*Tdiff';
end

end