# Dynamix

> Chromatin Stripe Calling Pipeline using a Bayesian NDLM

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

## Output

Once execution is complete, you can find the chromatin stripes in `MyDir/chr12_stripes.bedpe`.

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
   - Larger values yield a smoother NDLM smoothing mean distribution of the 1D signal probability sum (denoted as `sm` in Soomin Han).
5. **var_prior_ndlm** `1e4`
   - The prior variance for the initial state of the signal probability sum in the Bayesian NDLM.
6. **threshold_rowsums_em_p** `0.5`
   - The threshold for the intensity probability of the 1D signal probability sums to determine stripe anchor candidate genomic indices.
7. **peak_prominence** `0.1`
   - The peak prominence parameter for finding local maxima in the NDLM smoothing distribution of 1D signal probability sums.
8. **threshold_mean_to_peak** `0.1`
   - The ratio of (probability mass from the NDLM smoothing mean to the cut value) to (probability mass from the NDLM smoothing mean to the 1D signal probability sum).
   - The `cut_value` is used to determine the stripe width (denoted as `ε` in Soomin Han).
9. **max_stripe_width** `3e5`
   - The maximum stripe width to call.
10. **threshold_baseline_p** `0.95`
    - The baseline score to subtract from the score vector to determine stripe length (denoted as `δ` in Soomin Han).
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
To modify other parameters, open and update the values in `My_R_Library_Path/Dynamix/params.txt`, then save the file.
