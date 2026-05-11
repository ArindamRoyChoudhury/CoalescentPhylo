# Package ‘CoalescentPhylo’

## Overview

`CoalescentPhylo` is an R package for phylogenetic tree estimation through the method of root distances from allele count data using the coalescent model.
## Installation

Install the released version of this package from CRAN

```r
install.packages("CoalescentPhylo")
```

Or install the development version of this package from GitHub

```r
remotes::install_github('ArindamRoyChoudhury/CoalescentPhylo', upgrade="never")
```


## Example

```r
library(CoalescentPhylo)

# Load built-in example dataset (9 taxa x 2000 loci)
data(Human_Allele_Count_Data)

# Inspect dimensions: rows = taxa, columns = loci
dim(Human_Allele_Count_Data)

# Preview first few rows and columns
Human_Allele_Count_Data[1:3, 1:5]

# Estimate phylogenetic tree using the Coalescent Model
tree <- RD(mat_allele_count = Human_Allele_Count_Data, n.cores = NULL)

# Summarize the result
print(tree)

# Plot the labeled phylogenetic tree
plot(tree$labeled_tree)
```
