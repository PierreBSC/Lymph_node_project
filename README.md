# Lymph Node project
List of all R codes used for the analysis of single cell data in the paper Blecher et al.
Data have not been uploaded yet to GEO but they will be soon.

Each 'analysis' script corresponds to a specific part of the paper :

- Script_total_cells_analysis.R : used to perform the analysis of the total cell fraction. Corresponds to Figure 1 and 3.
- Script_Antigen_pos_analysis.R : used to perform the analysis of the antigen positive cell fraction. Corresponds to Figure 2.
- Script_all_mix_cells_analysis.R : used to perform the analysis of the merged total cell and antigen positive fractions. Corresponds to Figure 2 and 3.
- Script_IFN_NK_depletion.R : used to perform the analysis of the NK/IFNg depletion experiment. Corresponds to Figure 4.


'Supplementary analysis' scripts corresponds to analysis not shown in the main figures of the paper or analysis asked by the reviewers :

- Supplementary_analysis_total_cells.R : additional QC for total cell fraction.
- Supplementary_analysis_ag_pos_cells.R : additional QC for antigen positive cell fraction.
- Supplementary_analysis_all_mix_cells.R : detailed study of migratory dendritic cells diversity.

Other scripts were used to preprocess the Index Sorting Data :

- AFA_pipeline_script.R : script containing several functions based on FlowCore package to extract, clean and map index sorting data.
- Index_sorting_total_cells.R : pre-processing of Index Sorting data for total cell fraction.
- Index_sorting_Antigen_pos.R : pre-processing of Index Sorting data for antigen positive cell fraction.
