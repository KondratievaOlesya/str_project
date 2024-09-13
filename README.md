# Expression Analysis across different Colorectal Tumors

## Introduction

**Short Tandem Repeats (STRs)** are small sections of DNA where the same sequence of nucleotides (like A, T, C, G) repeats over and over again. For example, **"ATATATAT"** is a STR with **"AT"** motif repeated 4 times. 

STRs are tricky to analyze because they are repetitive and can easily change, which makes them difficult to track with standard tools. However, recent research shows that some STRs can actually influence how much a gene is expressed. These are called **expression STRs (eSTRs)**. 

We are focusing on colorectal cancer, a type of cancer that can be divided into groups based on **microsatellite instability (MSI)**. MSI happens when the mismatch repair system in DNA is deficient, leading to lots of mutations in DNA sequences, including STRs. 


### Literature
1. Presentation with characteristics of MSI-H, MSS/TMB-H and MSS/TMB-L [[link](https://docs.google.com/presentation/d/1btpGaDmqaFcaVvUyUpp_v7kEHXu7n5NfMNS083s5adA/edit?usp=sharing)]
2. eSTRs in colorectal cancer. [[paper](https://doi.org/10.1038/s41598-024-53739-0)]

    Verbiest, M.A., Lundström, O., Xia, F. et al. Short tandem repeat mutations regulate gene expression in colorectal cancer. Sci Rep 14, 3331 (2024). https://doi.org/10.1038/s41598-024-53739-0




## Key questions:

1. Which eSTRs are associated with MSS/TMB-H? 

2. Are there specific STRs that are hypermutated only in MSS/TMB-H?

3. What is the potential of MSS/TMB-H to respond to immunotherapy based on STR associations?

## Workflow

#### Differential Expression (DE) analysis for gene and eSTR expression

**Objective:** Differential expression analysis to compare gene expression and eSTR expression between MSS/TMB-H and other subtypes (MSS/TMB-L and MSI).

**Questions to address:**

1. Are there eSTRs that assosiated with MSS/TMB-H or other groups?

2. What are the role of genes that assosiated with MSS/TMB-H? Is there anything related to immune activity?

**Tools:**
R (preferred): Packages like DESeq2 and EdgeR.


#### STR mutation differences across groups

**Objective:** Determine how STR mutations differ across MSS/TMB-H, MSS/TMB-L, and MSI, building on the known differences between MSS and MSI.

**Questions to Address:**

1. Is there a clear mutation pattern in STRs that distinguishes MSS/TMB-H from MSS/TMB-L tumors?

**Tools:**
R or Python for mutation rate comparisons and clustering analyses.

