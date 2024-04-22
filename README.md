# ckd-workflow
Repository for "Variants in tubule epithelial regulatory elements mediate most heritable differences in human kidney function" paper.

If you'd like to run the notebooks, install the included libraries, and update all the file paths in them.

The workflow requires Seurat_4.3. Other library versions are listed in R_session_info.txt

## Primary analysis
* 00_cellranger_count 
Runs Cellranger count on the fastq files. 

Input
- S1_resources/L?_library_reference.csv 
- Cellranger transcriptome reference, refdata-gex-GRCh38-2020-A
- Guide reference: S1_resources/guide_reference.csv

Output
- Cellranger reports 
- Feature-barcode matrix
- Guide calls (not used)

* 01_souporcell_demux.sh
Labels each cell by the donor. 

Input 
- barcodes.tsv.gz from Cellranger
- Genome reference 

Output
- Cells labeled by donor IDs


## Secondary analysis
* **00_QC.R.ipynb**: 
Single cell QC to eliminate doublets, low-count cells, and to inspect the data.

* **01_Integrate.R.ipynb**:
Integrates data from multiple 10x lanes into a single Seurat file.

* **02_Post-integration_QC_and_viz.ipynb**:
Visualization of the integrated data for inspection.

* **03_Guide_calling.ipynb**:
Guide calling using the Poisson-Gaussion model from *Scalable single-cell CRISPR screens by direct guide RNA capture and targeted library enrichment", Replogle JM et al., Nature Biotechnology 2020*

* **03b_Post-guide_calling_QC_and_viz.ipynb**:
Visual inspection of guide calls.

* **04_DE_by_vector.ipynb**:
Differential expression testing separately for each dual-guide vector using single cell methods.

* **04_DE_by_vector_pseudobulk.ipynb**:
Differential expression testing separately for each dual-guide vector using pseudobulk methods.

* **04_DE_by_target_pseudobulk.ipynb**:
Differential expression testing separately by combining dual-guide vectors that target the same gene/DE, and using pseudobulk methods.
 
* **04_DE_by_target.ipynb**:
Differential expression using testing single cell methods by combining the vectors associated with each target gene or distal element.

* **04_DE_by_target_genomewide.ipynb**:
Differential expression using testing single cell methods by combining the vectors associated with each target gene or distal element. Tests all genes. (Other notebooks test only genes within 1Mb up/downstream of the target.)

* **04_DE_by_donor.ipynb**:
Differential expression testing using single cell methods run separately on each donor's cells.

