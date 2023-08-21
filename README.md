# Cross-Population-Portability

MC-ANOVA is an extension of HD-ANOVA [Vazquez et al., 2020](https://pubmed.ncbi.nlm.nih.gov/33315963/) that predicts the portablity of SNP segments in the context of cross-population Polygenic Risk Scores (PGS). The goal is to estimate the extent of genome differentiation with within and across population R-squared. The [MC_ANOVA.R](https://github.com/lupiA/Cross-Population-Portability/blob/main/MC-ANOVA.R) function draws genetic values for QTL from a local core of SNPs and then predicts those values using SNPs not in the core or randomly chosen to be QTL. The R-squared is the squared correlation between the generated genetic values and the predicted values.
\
\
We also provide a function, [getSegments.R](https://github.com/lupiA/Cross-Population-Portability/blob/main/getSegments.R), to group SNPs into local segments based on a provided Kbp size (e.g., 10 Kbp) and minimum number of SNPs (e.g., 10 SNPs).
\
\
Finally, we have provided an interactive tool, an [R Shiny App](https://github.com/lupiA/Cross-Population-Portability/blob/main/R-shiny-app), in which users can input a single SNP (base pair [BP] position), range of SNPs (BP positions), or a comma-separated list of SNPs, and the App will output portability and marker information. Users are able to download the main output from the App to a .csv file.

## Portability Pipeline
1. MC-ANOVA
2. Portability Map
3. GWAS
4. Marker effects (Bayes)
5. PGS

### Example
1. MC-ANOVA
\
After downloading the [MC_ANOVA.R](https://github.com/lupiA/Cross-Population-Portability/blob/main/MC-ANOVA.R) and [getSegments.R](https://github.com/lupiA/Cross-Population-Portability/blob/main/getSegments.R) functions and the [geno_map_example.csv](https://github.com/lupiA/Cross-Population-Portability/blob/main/geno_map_example.csv) file (requires the R package [BGData](https://github.com/QuantGen/BGData/tree/master)):

```
# Load necessary packages
library(BGData)

# Set seed for reproducibility
set.seed(12345)

# Generate genotypes (100 subjects and 500 SNPs)
n <- 100
p <- 500
genotype_map <- read.csv("~/geno_map_example.csv", header = TRUE)
X <- matrix(sample(0:2, n * p, replace = TRUE), ncol = p)
colnames(X) <- genotype_map$SNPs

# Assign population IDs (80% to population 1, 20% to population 2)
n_1 <- 0.8 * n
n_2 <- 0.2 * n
population <- sample(c(rep("Pop_1", round(n_1)), rep("Pop_2", round(n_2))))
rownames(X) <- population

# Initialize MAP and define 10 Kbp segments (>9 SNPs per segment)
minSNPs <- 10
minBP <- 10e3
MAP <- genotype_map
MAP$segments <- getSegments(MAP$base_pair_position, chr = MAP$chromosome, minBPSize = minBP, minSize = minSNPs, verbose = TRUE)

# Perform MC-ANOVA
# Running one central segment with 10 SNP flanking buffer
core <- which(MAP$segments == 5)
flank_size <- 10
chunk_start <- min(core) - flank_size
chunk_end <- max(core) + flank_size
chunk <- chunk_start:chunk_end
isCore <- chunk %in% core

X_1 <- X[rownames(X) == "Pop_1", chunk]
X_2 <- X[rownames(X) == "Pop_2", chunk]

# Set parameters for MC-ANOVA
lambda <- 1e-8
nRep <- 300
nQTL <- 3

# Perform MC-ANOVA analysis
out <- MC_ANOVA(X_1, X2 = X_2, core = which(isCore), lambda = lambda, nQTL = nQTL, nRep = nRep)

# Extract R-squared estimates
R_squared_within <- out[1, 1]
R_squared_across <- out[2, 1]
```
