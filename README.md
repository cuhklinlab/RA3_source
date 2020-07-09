# The source code of RA3
The source code for the reproduction of results in "A reference-guided approach for epigenetic characterization of single cells". <br/>
We also provide a [R-based implementation](https://github.com/cuhklinlab/RA3) for the convenience of the community.

# Requirements
- Matlab R2019a for the RA3 model
- Python 3.6.10 for cell clustering
- R 3.6.1 for data visualization and downstream analyses

# Usage instructions
|Code|Function|
|:-|:-|
|`RA3_func.m`|main function of RA3 model| 
|`RA3_script.m`|data preprocessing and model initialization of RA3|
|`run_RA3_script.m`|running script of RA3|
|`plot_case_fig1.R`|Fig. 1bcef|
|`plot_residual.R`|Fig. 1d|
|`plot_case.R`|Fig. 2abc, Supplementary Fig. 1bcde, and Supplementary Fig. 2abdef|
|`plot_case_pbmc.R`|Fig. 2d and Supplementary Fig. 2cg|
|`plot_fig3_scatter.R`|Fig. 3|
|`run_clustering_RA3.py`|clustering of RA3|
|`plot_fig3-2_AMI_bar.R`|Fig. 4a and Supplementary Fig. 4a|
|`plot_fig4_AMI_line.R`|Fig. 4b and Supplementary Fig. 4b|
|`plot_bar_ragi.R`|Fig. 4c|
|`trajectory.R`|trajectory inference, Fig. 5a|
|`cluster_specific_peak.R`|cluster specific peaks|
|`chromVAR_motif.R`|motif analysis, Fig. 5b|
|`plot_corr3.m`|Supplementary Fig. 1a|
|`plot_umap.R`|Supplementary Fig. 3|
|`plot_GmHekInSiMCA.R`|Supplementary Fig. 5|
|`plot_dropout_cell_sum.R`|Supplementary Fig. 6|
|`plot_case_pbmc_all_other.R`|Supplementary Fig. 7|
|`plot_supp_diff_dim.R`|Supplementary Fig. 8|

# License
This project is built under license **GNU GENERAL PUBLIC LICENSE (GPL)**.