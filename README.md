# Cross-Population-Portability

MC-ANOVA is an extension of HD-ANOVA [Vazquez et al., 2020](https://pubmed.ncbi.nlm.nih.gov/33315963/) that predicts the portablity of SNP segments in the context of cross-population Polygenic Risk Scores (PGS). The goal is to estimate the extent of genome differentiation with within and across population R-squared. The MC-ANOVA.R (LINK HERE) function draws genetic values for QTL from a local core of SNPs and then predicts those values using SNPs not in the core or randomly chosen to be QTL. The R-squared is the squared correlation between the generated genetic values and the predicted values.
\
\
We also provide a function, getSegments.R (LINK HERE), to group SNPs into local segments based on a provided Kbp size (e.g., 10 Kbp) and minimum number of SNPs (e.g., 10 SNPs).
\
\
Finally, we have provided an interactive tool, an R Shiny App (LINK HERE), in which users can input a single SNP (base pair [BP] position), range of SNPs (BP positions), or a comma-separated list of SNPs, and the App will output portability and marker information. Users are able to download the main output from the App to a .csv file.

\
## Portability Pipeline

1. MC-ANOVA
2. Portability Map
3. GWAS
4. Marker effects (Bayes)
5. PGS

### Example

#### 1. MC-ANOVA
After downloading the MC-ANOVA.R function (LINK HERE):

```
set.seed(12345)

# generate genotypes (100 subjects and 500 SNPs)
n <- 100
p <- 500
genotype_data <- matrix(sample(0:2, n * p, replace = TRUE), ncol = p)

# assign population IDs (80% to population 1, 20% to population 2)
n_1 <- 0.8 * n
n_2 <- 0.2 * n
population <- sample(c(rep("Pop_1", round(n_1)), rep("Pop_2", round(n_2))))


```
