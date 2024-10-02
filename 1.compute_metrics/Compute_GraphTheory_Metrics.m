%% COMPUTE GRAPH THEORY METRICS
% This script computes various graph theory metrics on subjects' 
% structural and functional connectivity matrices and saves the metrics 
% to the specified output directory. The computed metrics include:
% 
% 1. Community Structure and Modularity
% 2. Local Clustering Coefficient
% 3. Weighted Degree
% 4. SC-FC Coupling Metric

% Before running this script, ensure that you have:
% 1. Specified the paths for the repository, data, and output directories.
% 2. The `community_louvain` and `clustering_coef_wu` functions available 
%    in the MATLAB path. They should be located in the 'Rokos2024_SCFC_NetworkAnalyses/1.compute_metrics/functions' directory.

%% Step 1: Specify Directories & Paths

% Define base directories
repo_dir = '/path/to/directory/Rokos2024_SCFC_NetworkAnalyses/'; % MODIFY path to the cloned repository
data_dir = fullfile(repo_dir, 'data/'); % Data directory
outputs_dir = fullfile(repo_dir, 'outputs/'); % Graph theory metric directory

% Add repository directory to MATLAB path
addpath(genpath(repo_dir));

%% Step 2: Load in All Matrices
load([data_dir, 'subject_matrices.mat'])

%% Step 3: Compute the Graph Theory Metrics by running each cell
%% SC Community
SC_community=[];
SC_M=[];
gamma = 1;
for s=1:78
    W = SC_matrices{s};

    %Calculate
    [M, Q] = community_louvain(W,gamma);
    
     %Append Q metric  
     SC_community =[SC_community; Q];

    %Append partition
     SC_M=[SC_M; M'];
end

%% FC Community
FC_community=[];
FC_M=[];
gamma = 3;
for s=1:78
    %Set Diagonal to 0
    A=FC_matrices{s};
    n=size(A,1);
    A(1:n+1:n*n)=0;
    
    %Fisher-Z Transform
    matrix = atanh(A);
    
    %Calculate
    [M, Q] = community_louvain(matrix,gamma,'','negative_asym');
    
     %Append Q metric     
     FC_community =[FC_community; Q];
    
    %Append partition
     FC_M=[FC_M; M'];
end

%% SC Local Clustering
SC_clustering=[];

for s=1:78
    W = SC_matrices{s};

    %Calculate
    C = clustering_coef_wu(W);
    
    %Append LC metric  
    SC_clustering(s,:) = C;

end

%% FC Local Clustering
FC_clustering=[];

for s=1:78
        %Set Diagonal to 0
        A=FC_matrices{s};
        n=size(A,1);
        A(1:n+1:n*n)=0;
        
        %Fisher-Z Transform
        W = atanh(A);
        
        %Calculate
        [C_pos,C_neg,Ctot_pos,Ctot_neg] = clustering_coef_wu_sign(W,3);
        
        %Append LC metric  
        FC_clustering(s,:) = C_pos;   
end

%% Weighted Degree & SC-FC Coupling
    SCFC_coupling = [];
    SC_weighted_degree = [];
    FC_weighted_degree = [];
    
    for s = 1:78
        %Get Input Matrices:
        subSC = SC_matrices{s};
        subFC = FC_matrices{s};
        
        %Set Diagonal to 0
        A=FC_matrices{s};
        n=size(A,1);
        A(1:n+1:n*n)=0;
        
        %Fisher-Z Transform
        matrix=A;
        FC_z = atanh(matrix);

        %Calculate Column Wise Average Vectors
        subSC_col = mean(subSC, 2);
        subFC_col = mean(FC_z, 2);

        %Create SC Stacked Weighted Degrees
        SC_weighted_degree = [SC_weighted_degree; subSC_col'];

        %Create FC Stacked Weighted Degrees
        FC_weighted_degree = [FC_weighted_degree; subFC_col'];

        %Calculate correlation of SC & FC weighted degrees:
        tmp_corr = zeros(1, 193); % Initialize tmp_corr vector
     
        for region = 1:193
            tmp_corr(region) = corr(subSC(:,region), FC_z(:,region), 'type', 'Spearman');
        end
        
        SCFC_coupling = [SCFC_coupling; tmp_corr]; % Append correlation value of SC & FC

    end

%% Step 4: Save Metrics to Output Directory
save(fullfile(outputs_dir, 'graph_theory_metrics.mat'), ...
    'SC_community', 'SC_M', 'FC_community', 'FC_M', ...
    'SC_clustering', 'FC_clustering', ...
    'SC_weighted_degree', 'FC_weighted_degree', 'SCFC_coupling');
