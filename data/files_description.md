# Data Sources
## Clinical information

1. **coad_clinical.csv:** Clinical information for colorectal tumors from the TCGA-COAD dataset.
2. **read_clinical.csv:** Clinical information for rectal tumors from the TCGA-READ dataset.

These files contain clinical data that we will use to categorize tumors for further analysis.

### Description of CSV Files

The two clinical files share similar structures and contain the following important columns:

1. Tumor_Sample_Barcode: A unique identifier for each tumor sample.
2. MSI_TMB: A combined column that contains both the MSI status (Microsatellite Instability) and TMB value (Tumor Mutational Burden). This column is categorizing tumors into the following three groups:
    * MSI-H     (microsatellite instability-high)
    * MSS/TMB-H (microsatellite stable with high tumor mutational burden)
    * MSS/TMB-L (microsatellite stable with low tumor mutational burden)

## Gene expression
You can download it from [google drive](https://drive.google.com/file/d/1i3CHfvPARuu9rwnucBT3KimP9p-bNUxc/view?usp=sharing)

In the `gene_expression` folder, we have gene expression data for both TCGA-COAD and TCGA-READ cohorts, each located in their respective subfolders:

1. `coad/`: Contains gene expression data for colorectal tumors.
2. `read/:` Contains gene expression data for rectal tumors.

Each subfolder includes several files with different types of normalization:
1. `counts.csv`: Contains the raw gene expression counts for each sample.
2. `fpkm.csv:` Contains gene expression data normalized using FPKM (Fragments Per Kilobase of transcript per Million mapped reads).
3. `fpkm_uq.csv:` Contains Upper Quartile (UQ) FPKM normalized data, which adjusts for sequencing depth and variability across samples.
4. `tmm.csv:` Contains gene expression data normalized using the TMM (Trimmed Mean of M-values) method, often used to correct for library size differences.
5. `tpm.csv:` Contains gene expression data normalized by TPM (Transcripts Per Million), a method that accounts for gene length and sequencing depth.