% Table 1 (Discrete: Tauchen)
start1 = [0, 10, 50, 90, 95, 99];
mean1 = [0.14679, 0.041005, -0.04104, -0.10807, -0.14588, -0.34201];
std1 = [0.2686, 0.2323, 0.23271, 0.23845, 0.24549, 0.38598];

% Table 2 (Discrete: Rouwenhorst)
start2 = [0, 10, 50, 90, 95, 99];
mean2 = [0.13549, 0.040995, -0.040534, -0.10224, -0.14161, -0.27319];
std2 = [0.26146, 0.22812, 0.22805, 0.23485, 0.25026, 0.37611];

% Table 3 (Data)
start3 = [0, 10, 50, 90, 95, 99];
mean3 = [0.16, 0.02, -0.03, -0.08, -0.12, -0.23];
std3 = [0.41, 0.2, 0.21, 0.25, 0.3, 0.48];

% Table 4 (Continuous: Rouwenhorst)
start4 = [0, 10, 50, 90, 95, 99];
mean4 = [0.13146, 0.031615, -0.035691, -0.085678, -0.11747, -0.23438];
variance4 = [0.06800, 0.037749, 0.039568, 0.060507, 0.071352, 0.11751];
std4 = sqrt(variance4);

% Masaya's Table
start5 = [0, 10, 50, 90, 95, 99];
mean5 = [0.12146, 0.041615, -0.025691, -0.075678, -0.10747, -0.21438];
variance5 = [0.05800, 0.027749, 0.029568, 0.050507, 0.061352, 0.10751];
std5 = sqrt(variance5); % Compute std from variance

% Create tiled layout for better spacing
pdf_file_name = 'growth_rates_plots.pdf';
figure;
t = tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact'); % Two rows, one column

% Plot Mean Growth Rates
nexttile;
plot(start1, mean1, '-o', 'DisplayName', 'Discrete: Tauchen'); hold on;
plot(start2, mean2, '-s', 'DisplayName', 'Discrete: Rouwenhorst'); hold on;
plot(start3, mean3, '-^', 'DisplayName', 'Data'); hold on;
plot(start4, mean4, '-d', 'DisplayName', 'Continuous: Rouwenhorst'); hold on;
plot(start5, mean5, '-p', 'DisplayName', "Masaya's Table"); hold on;
xlabel('Percentile Start');
ylabel('Mean Growth');
title('Mean Growth Rates');
legend('Location', 'best');
grid on;

% Plot Std Growth Rates
nexttile;
plot(start1, std1, '-o', 'DisplayName', 'Discrete: Tauchen'); hold on;
plot(start2, std2, '-s', 'DisplayName', 'Discrete: Rouwenhorst'); hold on;
plot(start3, std3, '-^', 'DisplayName', 'Data'); hold on;
plot(start4, std4, '-d', 'DisplayName', 'Continuous: Rouwenhorst'); hold on;
plot(start5, std5, '-p', 'DisplayName', "Masaya's Table"); hold on;
xlabel('Percentile Start');
ylabel('Standard Deviation of Growth');
title('Standard Deviation of Growth Rates');
legend('Location', 'best');
grid on;

% Save the tiled layout to a single PDF file
set(gcf, 'PaperOrientation', 'portrait', 'PaperSize', [8.5 11], 'PaperPosition', [0 0 8.5 11]);
print('-dpdf', pdf_file_name);

disp(['Plots saved to ', pdf_file_name]);
