# HSCA Extended – Integration and Annotation Scripts

This folder contains scripts for the integration and annotation of the **HSCA Extended**.

## Script Overview

- **1_merge_and_prep_extended_datasets.R**  
  Merges 21 processed single-cell datasets **not included in the HSCA core** (datasets without PSU cells) into a unified Seurat object.  
  Also removes gene symbol artifacts from the original studies.

- **2_scarches_HSCA.ipynb**  
  Maps HSCA Extension datasets onto the HSCA Core using **scArches**.

- **3_label_transfer_HSCA.ipynb**  
  Transfers cell type labels from the HSCA Core to the HSCA Extension datasets via the **scArches** framework.

- **4_anndata_to_rds.ipynb**  
  Saves the relevant information from the HSCA Extended Anndata file after transfer learning, enabling export as an RDS file for further analysis in **R**.

- **5_reconstruct_hsca_extended.R**  
  Reads the output from `4_anndata_to_rds.ipynb` to reconstruct the HSCA Extended object and performs basic processing to generate UMAP embeddings and clustering.

- **6_annotate_hsca_extended.R**  
  Annotates the final HSCA Extended across five hierarchical cell type levels.  
  Uses a bottom-up approach after subclustering major cell types.

- **helpers.R**  
  Contains utility functions used across the main scripts.

- **scripts_subclustering/**  
  Directory containing additional scripts for subclustering specific cell populations.

## Workflow

1. **Merge and preprocess** extended single-cell datasets → unified Seurat object  
2. **Map and transfer labels** using scArches → create joint embedding and draft cell type annotation  
3. **Export relevant data** to RDS → reconstruct HSCA Extended in R  
4. **Hierarchical annotation** → subcluster and assign final cell identities

