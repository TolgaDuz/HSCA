# HSCA Core – Integration and Annotation Scripts

This folder contains the scripts used to perform the integration and annotation of the **HSCA core**.  

## Script Overview

- **1_merge_raw_datasets_hsca_core.R**  
  Merges 13 processed single-cell datasets containing **PSU cells** into a unified Seurat object.  

- **2_integration_analysis_scvi.R**  
  Performs data integration with **scVI** and generates a UMAP embedding.  

- **3_annotate_and_build_hsca_core.R**  
  Annotates the final HSCA core across five hierarchical cell type levels.  
  Annotation is carried out using a **bottom-up approach** after subclustering of the major cell types.  

- **helpers.R**  
  Contains utility functions used across the main scripts.  

- **scripts_subclustering/**  
  Directory containing scripts for subclustering.  

## Workflow

1. **Merging** of single-cell datasets → Seurat object  
2. **Integration & UMAP** with scVI  
3. **Hierarchical annotation** via subclustering  
