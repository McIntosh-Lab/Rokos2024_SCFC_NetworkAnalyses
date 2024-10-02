# Analysis code for: "Examining Relationships between Functional and Structural Brain Network Architecture, Age, and Attention Skills in Early Childhood"

This repository contains the code used for analyses in:

Rokos, L., Bray, S.L., Neudorf, J., Samson, A.D., Shen, K., & McIntosh, A.R. (2024). "Examining Relationships between Functional and Structural Brain Network Architecture, Age, and Attention skills in Early Childhood"

## Overall Workflow

1. **Clone the Repository**
   - Clone this repository to your local computer.

2. **PLS MATLAB Code**
   - Clone the PLS repository or download the PLS MATLAB code from [McIntosh Lab PLS GitHub](https://github.com/McIntosh-Lab/PLS) to your local computer.

For the scripts in: `~/Rokos2024_SCFC_NetworkAnalyses/1.compute_metrics`, `~/Rokos2024_SCFC_NetworkAnalyses/2.analyses_contrasts/` and `~/Rokos2024_SCFC_NetworkAnalyses/3.visualise_brainplots/`, you will need to modify paths, file names, and the variables of interest within the scripts first. "%MODIFY" in those scripts indicates where needs to be changed.

3. **Compute Graph Theory Metrics**
   - In MATLAB, run the `Compute_GraphTheory_Metrics.m` script located at `~/Rokos2024_SCFC_NetworkAnalyses/1.compute_metrics/` to compute and save the necessary graph theory metrics.

4. **Run Mean-Centered PLS Analyses**
   - In MATLAB, run the `Figure_2_Mean_Centered_PLS.m` script located at `~/Rokos2024_SCFC_NetworkAnalyses/2.analyses_contrasts/` to obtain results for Figure 2.

5. **Run Behavioural PLS (bPLS) Analyses**
   - In MATLAB, run one of the following scripts located at `~/Rokos2024_SCFC_NetworkAnalyses/2.analyses_contrasts/`:
     - For Figure 3: Run `Figure_3_Behavioural_PLS.m`
     - For Figure 4: Run the `Figure_4A_Behavioural_PLS.m` and `Figure_4B_Behavioural_PLS.m` scripts
     - For Figure 5: Run `Figure_5_Behavioural_PLS.m`
     - For Figure 6: Run `Figure_6_SCFC_Coupling_Combined_bPLS.m`
     - Note: These scripts include code to run the behavioural PLS analyses and to plot the contrasts.

6. **Visualize Results in R**
   - In R, run the `Brain_Plots.r` script located at `~/Rokos2024_SCFC_NetworkAnalyses/3.visualise_brainplots/` to plot the bootstrap ratios on a brain.
   - Note: The exact BSRs may be slightly different than in the paper because of resampling

## Necessary Tools

- **MATLAB PLS Code:** [McIntosh Lab PLS GitHub](https://github.com/McIntosh-Lab/PLS)
- **R Packages:** `ggseg` and `ggsegSchaefer`
- **Brain Connectivity Toolbox:** [BCT Website](https://sites.google.com/site/bctnet) (Required functions are located in `~/Rokos2024_SCFC_NetworkAnalyses/1.compute_metrics/functions`)
