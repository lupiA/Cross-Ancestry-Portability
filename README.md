# Cross-Population-Portability

MC-ANOVA is an extension of HD-ANOVA [Vazquez et al., 2020](https://pubmed.ncbi.nlm.nih.gov/33315963/) that estimates the extent of genome differentiation with within and across population R-squared. MC-ANOVA draws genetic values for QTL from a local core of SNPs and then predicts those values using SNPs not in the core or randomly chosen to be QTL. The R-squared is the squared correlation between the generated genetic values and the predicted values.

### Pipeline

1. MC-ANOVA
2. MAP
3. GWAS
4. Bayes
5. PGS

### Example

#### 1. MC-ANOVA
After downloading the MC-ANOVA.R function (link to MC-ANOVA.R fn:

