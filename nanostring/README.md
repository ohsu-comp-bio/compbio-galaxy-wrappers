## Analysis of Nanostring nCounter Vantage Protein Whole Slide Assay
##### Level 1
Signals are measured using the MAX/FLEX nCounter® Analysis System (Manual/Info)
This is the data in /X/Histopathology/Histopathology\ Shared\ Resource/CLINICAL/Nanostring/output/<batch_name>/rawdata.txt

##### Level 2 data
Based on Level 1 data, batch correction is preformed using factor analysis with removal of the wanted variation (RUV) method (Molania et al., 2019). RUV is preformed using a combination of  the “negative” antibody controls, MmAb-IgG1 and MmAb-IgG2, the negative and positive Hybridization ERCC spikes and the total S6 and total Histone H3 antibody signals. Any variation associated observed within 20 replicates of cell line technical replicates will be removed. 
Batches correction is performed with a cohort of 57 metastatic breast cancers. Consisting of:
 - 31: Breast Cancer not on Treatment
 - 16: Triple Negative Breast Cancer not on Treatment
 - 6: Triple Negative Breast Cancer on Treatment
This is the data in /X/Histopathology/Histopathology\ Shared\ Resource/CLINICAL/Nanostring/output

##### Level 3 data
Based on Level 1 data, batch correction is performed as describe in Level 2 including only batches containing samples of interest are combined, removing batch associated variation present only within the specified batches. 
This is the data in /X/Histopathology/Histopathology\ Shared\ Resource/CLINICAL/Nanostring/output_combining

**When analyzing the RPPA data:** 
 - For a single sample profiled in a single batch or multiple samples on the a single batch, either Level 1 or Level 2 should be good for analysis.
 - For multiple samples profiled in multiple batches, Level 3 is ideal because
   batch effects could occur among different batches. Level 2 could be used
with caution with the assumption that batch effect between your samples is
represented in variation of the 20 technical replicates used in the validation
cohort. There is also the potential to unnecessarily remove variation assocated
with batches that are not of interest 


### Workflow

![Batch Processing Workflow](images/Canvas%201.png)

![Combining Batches/Samples Workflow](images/Canvas%202.png)

### Assay Methods
Assay Methods
Nanostring 3D Vantage Solid Tumor Panel is comprised of 27 antibodies, including 13 phosphorylated protein targets, specifically designed to interrogate the MAPK and PI3K/mTOR signaling pathways (Lee et al., 2018).  This multiplexed panel allows for the simultaneously quantification of multiple proteins from a single section of FFPE tissue. 

Four micron sections of formalin-fixed paraffin-embedded (FFPE) cancer cell lines (controls) or 18-21 gauge cores of solid tumors were fixed immediately (within 3 minutes) of the biopsy to ensure retention of phosphorylated markers. Cores were subjected to citrate-based antigen retrieval and incubated overnight with the cocktail of oligo-tagged antibodies. After washing, the oligo-tags were released by UV light (3 minutes on a UV lightbox) and quantified using the Nanostring nCounter system. A set of 6 FFPE cancer cell lines were selected as positive controls and included on every run to assess antibody and control performance and to correct for batch effects. Batch correction was preformed using Removal of Unwanted Variation (RUV) using the replicate positive controls to estimate the factors associated with batch effect (Molania et al., 2018). RUV parameters were optimized by measuring the consistency of replicate controls and careful evaluation of outliers to ensure validity of results. This assay was validated for clinical use in the Knight Diagnostics Laboratories. 

References
Lee, J., Geiss, G.K., Demirkan, G., Vellano, C.P., Filanoski, B., Lu, Y., Ju, Z., Yu, S., Guo, H., Bogatzki, L.Y., et al. (2018). Implementation of a Multiplex and Quantitative Proteomics Platform for Assessing Protein Lysates Using DNA-Barcoded Antibodies. Mol Cell Proteomics 17, 1245–1258.
Molania, R., Gagnon-Bartsch, J.A., Dobrovic, A., and Speed, T.P. (2018). A new normalization for the Nanostring nCounter gene expression assay. BioRxiv.
 
### Versions

nCounter Digital Analyzer 5s
GeneRLF: NS_FFPE_STPath_Protein_v1.0
SystemAPF: 'n6_vDV1'
Output RCC FileVersion: 1.7
Nanostring SoftwareVersion: 4.0.1.8

Pipeline is currently run locally
 - Platform: x86_64-apple-darwin 15.6.0 (64-bit)
 - R version 3.5.2 (2018-12-20)
 - Python version Python 3.7.7
