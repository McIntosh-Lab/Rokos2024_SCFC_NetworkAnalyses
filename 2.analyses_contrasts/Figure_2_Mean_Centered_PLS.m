%% Mean Centered PLS Analysis and Saving BSR Values
% This script performs a Mean Centered Partial Least Squares (PLS) analysis 
% and saves the Bootstrap Ratio (BSR) values for further plotting.

% Before running this script, ensure that you have completed the following:
% 1) Computed the necessary brain metrics (e.g., FC_clustering, FC_weighted_degree, SCFC_coupling) 
%    using the `Compute_GraphTheory_Metrics.m` script.
% 2) Downloaded the PLS code (https://github.com/McIntosh-Lab/PLS) and added it to the MATLAB path.
% 3) Modified the following parameters in this script:
%    Step 1: Specify Directories & Paths
%       a) Define the paths for the repository, data, and output directories.
%       b) Modify the file names based on the specific figure you are generating
%        (e.g., 'FIGURE2A_bsrs.csv', 'FIGURE2B_bsrs.csv', or 'FIGURE2C_bsrs.csv').

%    Step 2: Set the Brain Data
%       a) Set `datamat_lst{1}` to the brain metric of interest (e.g., FC_clustering, FC_weighted_degree, or SCFC_coupling).

%    Step 3: Extract and Save Significant BSR Values
%       a) Set `lv` to 2 if reproducing Figure 1C (leave as 1 otherwise).


%% Step 1: Specify Directories & Paths

% Add PLS code to MATLAB path
addpath(genpath('/path/to/directory/plscmd')); %MODIFY

% Define base directories
repo_dir = '/path/to/directory/Rokos2024_SCFC_NetworkAnalyses/'; % MODIFY path to the cloned repository
data_dir = fullfile(repo_dir, 'data/'); % Data directory
outputs_dir = fullfile(repo_dir, 'outputs/'); % Output directory

% File Names
MCPLS_results_file = 'FIGURE2A_MCPLS.mat'; % MODIFY file name for saving the PLS result based on the figure (2A, 2B, or 2C)
bsr_csv_filename = 'FIGURE2A_bsrs.csv'; % MODIFY file name for saving significant BSR values based on the figure (2A, 2B, or 2C)

% Add repository directory to MATLAB path
addpath(genpath(repo_dir));

%% Step 2: Set the Brain Metric of Interest & Run the PLS Analysis
load(fullfile(outputs_dir, 'graph_theory_metrics.mat')); %Load in the graph theory metrics
datamat_lst{1}=FC_clustering; % MODIFY this to FC_clustering, FC_weighted_degree, or SCFC_coupling

num_subj=[39];%vector with the sample size per group
num_cond=2; %number of conditions per group
option.method=1;
option.num_boot=500;%number of bootstrap resamples
option.num_perm=1000;%number of permutation samples
option.meancentering_type=0;%mean centering type

% Run the PLS Analysis
result_mc = pls_analysis(datamat_lst,num_subj,num_cond,option);
result_mc.perm_result.sprob % Display the p-values'

% Save analysis output
save(fullfile(outputs_dir, MCPLS_results_file), 'result_mc');


%% Step 3: Extract and Save Significant BSR Values
lv=1; % LV of interest (MODIFY to 2 if reproducing Figure 1C)
thresh = 2; %BSR threshold
bsr = result_mc.boot_result.compare_u(:,lv);
sig_bsr_idx = bsr > thresh | bsr < (thresh*-1); % Identify significant BSR values

% Create an array with significant BSR values, setting others to zero
sig_bsr = zeros(size(bsr));
sig_bsr(sig_bsr_idx) = bsr(sig_bsr_idx);

% Save the sig_bsr array to a CSV file for plotting in R
csvwrite(fullfile(outputs_dir, bsr_csv_filename), sig_bsr); 
