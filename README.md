# Dynamix

> Chromatin Stripe Calling Pipeline using a Bayesian Normaly Dynamic Linear Model (NDLM)
> 
> ![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24079246.svg)

Dynamix identifies chromatin stripes from Hi-C contact files by applying a Bayesian Normal Dynamic Linear Model (NDLM) to signal probability matrices. The signal probability matrices are constructed using a zero-inflated negative binomial and negative binomial (ZINB+NB) mixture expectation-maximization (EM) algorithm.

> Soomin Han. (2026). *Dynamix identifies accurate chromatin stripes by Bayesian Normal Dynamic Linear Model on Hi-C signal probabilities*.

---

## Installation

To install Dynamix from GitHub, you first need to install the required Bioconductor dependencies for extracting GC content and mappability scores from restriction enzyme cutsites. 

Please run the following code in your R console:

```R
# 1. Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 2. Install required Bioconductor dependencies and reference genome
BiocManager::install(c(
    "BSgenome.Hsapiens.UCSC.hg19", 
    "BiocGenerics", 
    "IRanges", 
    "GenomicRanges", 
    "GenomeInfoDb", 
    "Biostrings", 
    "BSgenome", 
    "rtracklayer"
))

# 3. Install Dynamix from GitHub
devtools::install_github("soominhan511/Dynamix")
```

---

## Data Preparation

1. **Hi-C Data:** Prepare your `.hic` file.
2. **Mappability File:** Download the mappability bigWig file (`wgEncodeCrgMapabilityAlign50mer.bigWig`) from the UCSC Genome Browser:
   * [UCSC Mappability Download Link](https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/)

---

## Quick Start (Example: Chromosome 12)

Here is a quick example to call chromatin stripes on chromosome 12:

```R
library(Dynamix)

# Execute stripe calling
Call_Stripes(
    hic_file = 'MyDir/my_sample.hic', 
    wg_file = 'MyDir/wgEncodeCrgMapabilityAlign50mer.bigWig',
    output_path = 'MyDir/', 
    chr = 'chr12'
)
```

---

## Output File

#### 1. Stripe File
- **File Name:** `MyDir/chr12_stripes.bedpe`
- **Description:** This file includes all pixel interactions (i,j) located inside identified chromatin stripes.

| chr | start1 | end1 | start2 | end2 | stripe_number | stripe_width_range | stripe_length_range | stripe_width | stripe_length | p_value_stripe | direction | count | mu | signal_p | p_value_NB_MIX | pass |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| chr12 | 520001 | 530000.0 | 660001 | 670000.0 | 1 | 520000-610000 | 660000-2000000 | 90000 | 1340000 | 7.10e-05 | 5_prime | 35 | 127.33 | 0.0849 | 0.8918 | True |
| chr12 | 520001 | 530000.0 | 670001 | 680000.0 | 1 | 520000-610000 | 660000-2000000 | 90000 | 1340000 | 7.10e-05 | 5_prime | 33 | 122.40 | 0.0830 | 0.8951 | True |
| chr12 | 520001 | 530000.0 | 680001 | 690000.0 | 1 | 520000-610000 | 660000-2000000 | 90000 | 1340000 | 7.10e-05 | 5_prime | 23 | 111.53 | 0.0596 | 0.9325 | True |
| chr12 | 520001 | 530000.0 | 690001 | 700000.0 | 1 | 520000-610000 | 660000-2000000 | 90000 | 1340000 | 7.10e-05 | 5_prime | 27 | 107.94 | 0.0759 | 0.9070 | True |
| chr12 | 520001 | 530000.0 | 700001 | 710000.0 | 1 | 520000-610000 | 660000-2000000 | 90000 | 1340000 | 7.10e-05 | 5_prime | 26 | 95.97 | 0.0837 | 0.8940 | True |
| chr12 | 520001 | 530000.0 | 710001 | 720000.0 | 1 | 520000-610000 | 660000-2000000 | 90000 | 1340000 | 7.10e-05 | 5_prime | 31 | 93.61 | 0.1055 | 0.8544 | True |
| chr12 | 520001 | 530000.0 | 720001 | 730000.0 | 1 | 520000-610000 | 660000-2000000 | 90000 | 1340000 | 7.10e-05 | 5_prime | 37 | 91.21 | 0.1317 | 0.8023 | True |
| chr12 | 520001 | 530000.0 | 730001 | 740000.0 | 1 | 520000-610000 | 660000-2000000 | 90000 | 1340000 | 7.10e-05 | 5_prime | 43 | 85.84 | 0.1638 | 0.7335 | True |
| chr12 | 520001 | 530000.0 | 740001 | 750000.0 | 1 | 520000-610000 | 660000-2000000 | 90000 | 1340000 | 7.10e-05 | 5_prime | 25 | 82.81 | 0.0951 | 0.8738 | True |

**Column Description:**
- `chr`: Chromosome name
- `start1`: i * resolution (bp)
- `end1`: (i+1) * resolution (bp)
- `start2`: j * resolution (bp)
- `end2`: (j+1) * resolution (bp)
- `stripe_number`: Unique identifier for a stripe (all interactions within the same stripe share this value).
- `stripe_width_range`: Stripe width range (anchor) (all interactions within the same stripe share this value).
- `stripe_length_range`: Stripe length range (all interactions within the same stripe share this value).
- `stripe_width`: Stripe width (all interactions within the same stripe share this value).
- `stripe_length`: Stripe length (all interactions within the same stripe share this value).
- `p_value_stripe`: P value for stripe (all interactions within a stripe share this value)
- `direction`: `5_prime` for 5' stripes, `3_prime` for 3' stripes. (all interactions within the same stripe share this value).
- `count`: Raw Hi-C contact count for interaction (i,j)
- `mu`: Expected count of interaction (i,j)
- `signal_p`: Signal probability of interaction (i,j)
- `p_value_NB_MIX`: P value derived from the ZINB-NB mixture model
- `pass`: `True` if the stripe meets the filtering threshold (all interactions within a stripe share this value).

---

#### 2. Stat Files
- **2-1. File Name:** `MyDir/chr12_stats/chr12_pair_features_stats.csv`
  - **Information:** Contains the signal probabilities for all pairwise interactions (i,j) across chromosome 12.
- **2-2. File Name:** `MyDir/chr12_stats/ZINB_NB_mix_EM_coefficients.csv`
  - **Information:** Contains converged parameters from the ZINB and NB mixture models, including coefficient estimates, associated p-values, and AIC values.
- **2-3. File Name:** `MyDir/chr12_stats/chr12_rowsum_smoothed_ndlm_5_prime.png`
  - **Information:** Visualizes the Bayesian Normal Dynamic Linear Model (NDLM) smoothing means of 1D signal probability sums across genomic indices.
---

## Parameter Information and Default Values

If you execute the script as shown above, the stripes are obtained automatically with the following default parameters:

1. **resolution** `10000`
2. **window_size** `51e4`
   - This parameter is used to apply a direction filter to call 5' and 3' stripes.
   - To assign a direction filter (1, 2, or 3) to every interaction `(i, j)`, Dynamix compares the row-wise signal probability sum of the interaction `(i, j - window_size/resolution : j + window_size/resolution)` and the column-wise signal probability sum of the interaction `(i - window_size/resolution : i + window_size/resolution, j)`.
3. **down_sample_ratio** `0.01`
   - The ratio of the selected sample size to the total sample size.
   - Selected samples are included to converge ZINB+NB mixture EM model.
4. **smoothing_factor** `150`
   - The ratio of observational variance to system variance in the Bayesian NDLM.
   - Larger values yield a smoother NDLM smoothing mean distribution of the 1D signal probability sum (denoted as `sm` in [Han, 2026]).
5. **var_prior_ndlm** `1e4`
   - The prior variance for the initial state of the signal probability sum in the Bayesian NDLM.
6. **threshold_rowsums_em_p** `0.5`
   - The threshold for the intensity probability of the 1D signal probability sums to determine stripe anchor candidate genomic indices.
7. **peak_prominence** `0.1`
   - The peak prominence parameter for finding local maxima in the NDLM smoothing distribution of 1D signal probability sums.
8. **threshold_mean_to_peak** `0.1`
   - The ratio of (probability mass from the NDLM smoothing mean to the cut value) to (probability mass from the NDLM smoothing mean to the 1D signal probability sum).
   - The `cut_value` is used to determine the stripe width (denoted as `ε` in [Han, 2026]).
9. **max_stripe_width** `3e5`
   - The maximum stripe width to call.
10. **threshold_baseline_p** `0.95`
    - The baseline score to subtract from the score vector to determine stripe length (denoted as `δ` in [Han, 2026]).
11. **threshold_p_value_stripe** `5e-3`
    - The threshold p-value to filter passed stripes.
    - The `pass` column (e.g., in `chr12_stripes.bedpe`) becomes `TRUE` if the following conditions are met:
      - `p < (threshold_p_value_stripe)^2` for stripes of length shorter than 60 kbp.
      - `p < threshold_p_value_stripe` for stripes of length 60 kbp or longer.

---

## How to Modify Parameters

### 1. Modifying Resolution
To modify the resolution, include the `resolution` argument when calling the function (e.g., for `resolution = 5000`):

```R
Call_Stripes(
    hic_file = 'MyDir/my_sample.hic', 
    wg_file = 'MyDir/wgEncodeCrgMapabilityAlign50mer.bigWig',
    output_path = 'MyDir/', 
    chr = 'chr12', 
    resolution = 5000
)
```

### 2. Modifying Other Parameters (Except Resolution)
To modify parameters other than resolution, open and edit the configuration file located at your package installation path:

- `My_R_Library_Path/Dynamix/params.txt`

Modify the desired values and save the file before running the analysis.
