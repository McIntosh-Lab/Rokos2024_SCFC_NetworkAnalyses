%% SC-FC Coupling PLS Analysis and Plotting
% This script performs a behavioural Partial Least Squares (bPLS) analysis on the
% SCFC-coupling, SC Weighted Degree, and FC Weighted Degree brain data with age as the behavioral measure.
% The script then plots the contrast from the analysis.

%Before running this script, ensure that you have completed the following:
% 1) Computed the necessary brain metrics (i.e., SC_weighted_degree, FC_weighted_degree,
%    SCFC_coupling) using the `Compute_GraphTheory_Metrics.m` script.
% 2) Downloaded the PLS code (https://github.com/McIntosh-Lab/PLS) and added it to the MATLAB path.
% 3) Specified the paths for the repository, data, and output directories.
% 4) Modified the following parameters as needed:
%    a) `lv`: Set this to the significant LV in Step 3

%% Step 1: Specify Directories & Paths
% Add PLS code to MATLAB path
addpath(genpath('/path/to/directory/plscmd')); % MODIFY

% Define base directories
repo_dir = '/path/to/directory/Rokos2024_SCFC_NetworkAnalyses'; % MODIFY path to the cloned repository
data_dir = fullfile(repo_dir, 'data/'); % Data directory
outputs_dir = fullfile(repo_dir, 'outputs/'); % Output directory
figure_dir = fullfile(repo_dir, 'outputs/figures/'); % Figures directory

% File Names
bPLS_results_file = 'FIGURE6_bPLS.mat'; % MODIFY file name for saving bPLS results
figure_filename = 'FIGURE6.png'; % MODIFY file name for saving the contrast figure

% Add repository directory to MATLAB path
addpath(genpath(repo_dir));

%% Step 2: Define Brain and Behavioral Data
% Combine SC and FC metrics along with SCFC coupling for the analysis
load(fullfile(outputs_dir, 'graph_theory_metrics.mat')); %Load in the graph theory metrics
SCFC_datamat = [SC_weighted_degree; FC_weighted_degree; SCFC_coupling];

% Behavioral data: Age is used as the measure of interest
load(fullfile(data_dir, 'behaviouraldata.mat')); %Load in the behavioural data
behav = behaviouraldata(:,2); %
option.stacked_behavdata = [behav; behav; behav]; % Repeat for each brain metric

% Brain data matrix
datamat_lst{1} = SCFC_datamat;

% PLS analysis options
num_subj = [39];  % Sample size per group
num_cond = 6;     % Number of conditions per group (e.g., 3 metrics x 2 time points)
option.method = 3;                 % Behavioural PLS method
option.num_boot = 500;             % Number of bootstrap resamples
option.num_perm = 1000;            % Number of permutation samples
option.meancentering_type = 0;     % Mean centering type

% Run the PLS Analysis
result_SCFC = pls_analysis(datamat_lst, num_subj, num_cond, option);
disp(result_SCFC.perm_result.sprob); % Display p-values

% Save the result to a specified path
savePath = fullfile(outputs_dir, bPLS_results_file);
save(savePath, 'result_SCFC');

%% Step 3: Plotting SC-FC Coupling and Age Correlation
% This section plots the correlation structure between SC-FC coupling, SC Weighted Degree, FC Weighted Degree and age.

lv = 1;  % MODIFY LV of interest

% Define bar colors
barColors = [0 0.4470 0.7410; 0.6350 0.0780 0.1840; 0.4940 0.1840 0.5560];

% Create a figure for plotting
fig = figure('Units', 'inches', 'Position', [0, 0, 5.51, 5.51]);

% Extract data and plot with grouped colors
data = (result_SCFC.boot_result.orig_corr(:,lv)) * -1;
x = 1:6;
bar_handles = cell(3, 1);

for i = 1:numel(data)
    group_color = mod(floor((i - 1) / 2), 3) + 1;
    b = bar(x(i), data(i), 'FaceColor', barColors(group_color, :));
    bb(i) = b;
    hold on;
    bar_handles{group_color} = [bar_handles{group_color}, b];
end

% Define axis labels and title
axisLabelFontSize = 10;  % Font size for axis labels
titleFontSize = 12;      % Font size for the title

% Set X-axis labels and limits
set(gca, 'XTick', 1:size(result_SCFC.v(:,lv),1), 'XTickLabel', {'T1'; 'T2'; 'T1'; 'T2'; 'T1'; 'T2'}, ...
    'FontSize', axisLabelFontSize, 'LineWidth', 0.3);

% Define error bars
errhigh = ((result_SCFC.boot_result.ulcorr(:,lv) - result_SCFC.boot_result.orig_corr(:,lv)) * -1);
errlow = ((result_SCFC.boot_result.orig_corr(:,lv) - result_SCFC.boot_result.llcorr(:,lv)) * -1);

% Label the axes and title
xlabel('', 'FontWeight', 'bold', 'FontSize', axisLabelFontSize);
ylabel('Correlation between Age and Measures', 'FontWeight', 'bold', 'FontSize', axisLabelFontSize);
title(['LV=' num2str(lv) ', p=' num2str(result_SCFC.perm_result.sprob(lv))], 'FontSize', titleFontSize);

% Plot error bars
er = errorbar(x, data, errhigh, errlow);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 0.2;

% Save the figure
print(fullfile(figure_dir, figure_filename), '-dpng', '-r600');

hold off;
