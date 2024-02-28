# Gene-HWT

## Gene-based Hardy-Weinberg equilibrium test (gene-HWT) with genotype count data.

### Description
Perform gene-based Hardy-Weinberg equilibrium test (gene-HWT) using genotype count data and publicly available linkage disequilibrium (LD) information.

### Dependency

Gene-HWT functions depends on LDlinkR packages.

### Load and execute geneHWT.R
```r
source("geneHWT.R")
```

### Examples
● Here, using the `geneHWT` function, we perform gene-HWT for a gene (AK8 gene) using genotype count data for each SNPs on the gene. This function involves conducting a Hardy-Weinberg equilibrium test for each SNPs and integrating these results with the automatically obtained LD information (r^2) to obtain the z-value and corresponding P-value for gene-HWT.


View the test data (test_dt) generated by source("geneHWT.R").
```r
head(test_dt, 25)
    gene       rsID  AA Aa  aa
1    AK8  rs1000012  50 95  43
2    AK8 rs10117989  93 85  17
3    AK8 rs10901214  90 79  26
4    AK8 rs11243900  77 77  41
5    AK8 rs12000803   3 41 151
6    AK8  rs1410969 163 22   2
7    AK8  rs1609597 171 22   2
8    AK8 rs17149822 171 22   2
9    AK8   rs215156   6 39 150
10   AK8   rs215165  41 96  55
11   AK8   rs215192   7 46 141
12   AK8  rs2771994  10 62 123
13   AK8  rs2772007   4 48 143
14   AK8   rs415831  40 97  58
15   AK8  rs4962211  77 77  41
16   AK8  rs4962215   7 45 143
17   AK8  rs4962218  85 94  16
18   AK8   rs929293  92 85  18
19   AK8   rs932887   6 43 146
20   AK8   rs954390 174 20   1
21  NPNT rs10000111 142 46   5
22  NPNT rs11945032  99 83  13
23  NPNT rs13125117  62 95  38
24  NPNT rs17036413 129 60   6
25  NPNT rs17036425  15 15 111
```


Filter data for AK8 gene.
```r
test_dt_ak8 <- test_dt[test_dt$gene %in% c("AK8"),]
```

Perform gene-HWT using geneHWT function.
```r
geneHWT(dsin = test_dt_ak8, token = "******", pop = "EAS", genome_build = "grch37", n_r2 = 10)

```

Here, replace token = "**.." with your personal access token for accessing the LDlink API. For more details, refer to the LDlinkR documentation.



● Next, we demonstrate an example of performing gene-HWT using the geneHWT2 function, with both genotype count data for each SNPs in a gene and their corresponding linkage disequilibrium matrices prepared in advance.

Obtain linkage disequilibrium matrices (r^2) for variants (rsID) on AK8 gene using LDmatrix.

```r
ldmat <- LDmatrix(test_dt_ak8$rsID, pop = "EAS", r2d = "r2", token = "******", genome_build = "grch37")
ldmat <- ldmat[, -1]
```

Perform gene-HWT using geneHWT2 function.
```r
geneHWT2(dsin = test_dt_ak8, ldmat = ldmat)
```


## Functions

- **HWz Function:**
  - Arguments: AA, Aa, aa: Vectors indicating the number of individuals with each genotype (AA, Aa, aa).
  - Output: z: Vector indicating z-scores from HWT.



- **geneHWT Function:**
  - Arguments:
    - dsin: Data frame containing the following variables:
      - gene: Gene name
      - rsID: rsID of the gene variant
      - AA: Number of individuals with genotype AA
      - Aa: Number of individuals with genotype Aa
      - aa: Number of individuals with genotype aa
    - token: Personal access token for LDlinkR package. See LDlinkR package documentation for details.
    - pop: Population from the 1000 Genomes Project (e.g., YRI or CEU). See LDlinkR package documentation for details.
    - genome_build: Choose one of three options: for grch37 (hg19), use genome_build="grch37"; for grch38 (hg38), use genome_build="grch38"; for grch38 High Coverage (hg38) 1000 Genome Project dataset, use genome_build="grch38_high_coverage". Default is GRCh37 (hg19). See LDlinkR package documentation for details.
    - n_r2: Number of attempts to obtain r2 using LDmatrix. If attempts exceed n_r2, the function will abort.
  - Output: List containing the following variables:
    - gene: Gene name
    - n_vari: Number of variants (number of rsIDs) used in the analysis
    - HWz_com: z-values from geneHWT
    - HWp_com: P-values from geneHWT
 


- **geneHWT2 Function:**
  - Arguments:
    - dsin: Data frame containing the following variables:
      - gene: Gene name
      - rsID: rsID of the gene variant
      - AA: Number of individuals with genotype AA
      - Aa: Number of individuals with genotype Aa
      - aa: Number of individuals with genotype aa
    - ldmat: Matrix of linkage disequilibrium coefficients (r2) for rsIDs. Must be sorted according to the input rsIDs.
  - Output: List containing the following variables:
    - gene: Gene name
    - n_vari: Number of variants (number of rsIDs) used in the analysis
    - HWz_com: z-values from geneHWT
    - HWp_com: P-values from geneHWT
