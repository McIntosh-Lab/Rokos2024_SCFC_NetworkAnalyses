%% Behavioural PLS ANALYSES
% This script performs a Behavioural Partial Least Squares (bPLS) analysis 
% using brain data and behavioural measures, saves the significant 
% Bootstrap Ratio (BSR) values, and plots the bPLS constrasts for Figure 3. 
% This script also generates the scatter plot in Figure 3A. 

%Before running this script, ensure that you have completed the following:
% 1) Computed the necessary brain metric (i.e., SC_community) 
%   using the `Compute_GraphTheory_Metrics.m` script.
% 2) Downloaded the PLS code (https://github.com/McIntosh-Lab/PLS) and added it to the MATLAB path. 
% 3) Modified Script Parameters:
%    Step 1: Specify Directories & Paths
%      a) Define the paths for the repository, data, and output directories.
%      b) Modify the file names in the script based on the specific figure you are generating.

%   Step 2: Set the Brain and Behavioural Data
%      a) `datamat_lst{1}`: Ensure this is set to the brain metric of interest (i.e., SC_community).
%      b) `tasks`: Specify the behavioural data of interest. Options are:
%         1 = subject, 2 = age, 3 = selective_attention, 4 = sustained_attention, 5 = executive_attention.

%    Step 3: Adjust Plotting Details
%      a) `lv`: Ensure this is set to 1.
%      b) `ylim`: Adjust the y-axis limits based on the range of the data
%           (e.g., ylim([-1, 0]).
%      c) `xticklabels`: Modify this in Step 3 to reflect the specific behavioural measures you are plotting
%         (e.g., xticklabels = {'Age'; 'Sustained Attention'}).

%    Step 4: Extract and Save Significant BSR Values
%       a) Ensure LV of interest is set to 1.

%% Step 1: Specify Directories & Paths
% Add PLS code to MATLAB path
addpath(genpath('/path/to/directory/plscmd')); % MODIFY

% Define base directories
repo_dir = '/path/to/directory/Rokos2024_SCFC_NetworkAnalyses/'; % MODIFY path to the cloned repository
data_dir = fullfile(repo_dir, 'data/'); % Data directory
outputs_dir = fullfile(repo_dir, 'outputs/'); % Output directory
figure_dir = fullfile(repo_dir, 'outputs/figures/'); % Figures directory

% File Names
bPLS_results_file = 'FIGURE3B_bPLS.mat'; % MODIFY file name for saving bPLS results
figure_filename = 'FIGURE3B.png'; % MODIFY file name for saving the contrast figure
bsr_csv_filename = 'FIGURE3B_bsrs.csv'; % MODIFY file name for saving significant BSR values

% Add repository directory to MATLAB path
addpath(genpath(repo_dir));

%% Step 2: Set the Brain and Behavioural Data & Run the PLS Analysis
load(fullfile(outputs_dir, 'graph_theory_metrics.mat')); %Load in the graph theory metrics
datamat_lst{1} = SC_community; % MODIFY

% Specify the behavioural data of interest. Modify as needed.
load(fullfile(data_dir, 'behaviouraldata.mat')); %Load in the behavioural data
tasks = [2, 4]; % 1 = subject, 2 = age, 3 = selective_attention, 4 = sustained_attention, 5 = executive_attention
option.stacked_behavdata = behaviouraldata(:, tasks);

% PLS analysis options
num_subj = [39];  % Sample size per group
num_cond = 2;     % Number of conditions per group
option.method = 3;                  % Behavioural PLS method
option.num_boot = 500;              % Number of bootstrap resamples
option.num_perm = 1000;             % Number of permutation samples
option.meancentering_type = 0;      % Mean centering type

% Run the PLS Analysis
result_behavioural = pls_analysis(datamat_lst, num_subj, num_cond, option);
disp(result_behavioural.perm_result.sprob); % Display the p-values

% Save the result
save(fullfile(outputs_dir, bPLS_results_file), 'result_behavioural');


%% Step 3: Plotting Behavioural PLS Contrasts
% This section plots the correlation structure between brain data and 
% behavioural measures.

lv = 1; % MODIFY

% Create a figure for plotting
fig = figure('Units', 'inches', 'Position', [0, 0, 1.18, 1.18]);

% Define colors for the bars
barColors = [0.8 0.2 0.4; 0.4 0.7 0.9];

% Plot LV correlation structure
data = result_behavioural.boot_result.orig_corr(:, lv);
x = 1:numel(data);

% Create a bar plot with different colors for different timepoints
bar(x(1:2), data(1:2), 'FaceColor', barColors(1, :));
hold on;
bar(x(3:4), data(3:4), 'FaceColor', barColors(2, :));
ylim([-1, 0.5]); % Optional: MODIFY

% Add error bars
lower = result_behavioural.boot_result.orig_corr(:, lv) - result_behavioural.boot_result.llcorr(:, lv);
upper = result_behavioural.boot_result.ulcorr(:, lv) - result_behavioural.boot_result.orig_corr(:, lv);
er = errorbar(x, data, lower, upper);

er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 0.2;

% Customize axis labels and title
axisLabelFontSize = 3;
titleFontSize = 7;

xlabel('Behavioral Measures', 'FontWeight', 'bold', 'FontSize', axisLabelFontSize);
ylabel('Correlation with Behavioral Measures', 'FontWeight', 'bold', 'FontSize', axisLabelFontSize);
title(['LV=' num2str(lv) ', p=' num2str(result_behavioural.perm_result.sprob(lv))], 'FontSize', titleFontSize);

% Set X-axis labels
xticklabels = {'Age'; 'Sustained Attention'}; % MODIFY: X-axis labels
set(gca, 'XTick', 1:size(result_behavioural.v(:, lv), 1), 'XTickLabel', xticklabels, 'FontSize', axisLabelFontSize, 'LineWidth', 0.3);

% Save the figure with specified dimensions and resolution
print(fullfile(figure_dir, figure_filename), '-dpng', '-r600');

hold off;

%% Step 4: Extract and Save Significant BSR Values

lv = 1; % MODIFY
thresh = 2; % BSR threshold

% Extract BSR values for the specified LV
bsr = result_behavioural.boot_result.compare_u(:, lv);
sig_bsr_idx = bsr > thresh | bsr < -thresh;  % Identify significant BSR values

% Create an array with significant BSR values, setting others to zero
sig_bsr = zeros(size(bsr));
sig_bsr(sig_bsr_idx) = bsr(sig_bsr_idx);

% Save the sig_bsr array to a CSV file for plotting in R.
csvwrite(fullfile(outputs_dir, bsr_csv_filename), sig_bsr);

%% Step 5: Plot Figure 3A

figure;

% Plot scatter plot
scatter(behaviouraldata(:,2), SC_community, 'filled', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black');
hold on;

% Set y-axis limits
ylim([0.649, 0.69]);

% Calculate regression line
coefficients = polyfit(behaviouraldata(:,2), SC_community, 1);
x_range = [4.14,7.89];
y_range = polyval(coefficients, x_range);

% Plot the regression line
plot(x_range, y_range, 'r-', 'LineWidth', 2, 'Color', 'blue');

% Scale axis labels and title
axisLabelFontSize = 10;
titleFontSize = 12;
xlabel('Age', 'FontWeight', 'bold', 'FontSize', axisLabelFontSize);
ylabel('SC Modularity Statistic (Q)', 'FontWeight', 'bold', 'FontSize', axisLabelFontSize);
title('Age vs SC Modularity Statistic (Q)', 'FontSize', titleFontSize);

set(gca, 'FontSize', axisLabelFontSize, 'LineWidth', 0.3);

hold off;

print(fullfile(figure_dir, 'FIGURE3A'), '-dpng', '-r600');
hold off;

