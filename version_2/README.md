# MicroCellClust 2

Large-scale solver [2] implementing the MicroCellClust optimisation problem presented in [1]. 

### References

[1] A. Gerniers, O. Bricard and P. Dupont (2021). MicroCellClust: mining rare and highly specific subpopulations from single-cell expression data. *Bioinformatics*, 37(19), 3220-3227. https://doi.org/10.1093/bioinformatics/btab239

[2] Not yet published


## RScala

To run MicroCellClust in R, the user must install the `rscala` package (https://cran.r-project.org/web/packages/rscala/index.html), that can call the Scala implementation of MicroCellClust from within R. Once installed, he can check wether there is a version of Scala on his computer using:

``` R
rscala::scalaConfig()
```
This function tries to find Scala and Java on the user's computer and, if needed, downloads and installs Scala and Java in the user's `~/.rscala` directory.


### Initialising MicroCellClust

To run MicroCellClust, the user must first set his working directory (`setwd`) to the location of the `MicroCellClust`directory, and then load the functions contained in the `micro_cell_clust.R` script.
``` R
source("./src/main/R/micro_cell_clust.R")
```
Then, he must initialise an instance of RScala containing the MicroCellClust source code.
``` R
mcc.rs = initMccRscala()
```
This `mcc.rs` object needs to be passed as the first argument to any R function using Scala code.


## Preparing the data

MicroCellClust assumes the data is represented in a matrix (or data frame) containing positive values for presence of expression, and negative values in case of absence of (or negligible) expression. Such data can be obtained from (normalize) count data using an appropriate transformation, such as $\log_{10}(x + 0.1)$.

**Note:** Applying such a transformation means the data can no longer be stored in a sparse structure, which can be problematic when dealing with large structures of data. Therefore, all given `R` *will assume the given data are count values (or normalized counts)* and will apply the transformation internally. By default, they apply a $\log_{10}(x + 0.1)$ transformation. This behavior can be changed by setting the `data.tfo` parameter accordingly :

* `"log10"` : $\log_{10}(x + 0.1)$ transformation (default)
* `"log2"` : $\log_{2}(x + 0.5)$ transformation
* a custom `R` function with one parameter (corresponding to a single numerical value) that transforms the value passed as argument
* `"none"` if the data is already in the appropriate format and no transformation should be applied

Let's assume the (normalized) count data is contained in a matrix/dataframe called `mat` (assuming the rows represent genes and the columns cells). As discussed in the paper, MicroCellClust requires genes expressed in too many cells of the dataset to be removed.
``` R
gene.expr = rowSums(mat > 0) # Can be useful to plot a histogram to define a threshold. Typically keeping genes expressed in less that 25% keeps the vast majority of the genes
mat.filt = mat[gene.expr >= 2 & gene.expr <= 0.25 * ncol(mat), ]
```

MicroCellClust needs a vector containing the sum of positive expressions of each gene. This can be computed internally, but as it might be time-consuming for large instances, we recommend to pre-compute it once and store it in memory :
```R
gs = geneSum(mat.filt) # If needed, don't forget to set the `data.tfo` argument
```
**Remark:** By default, the R implementation represents the cells by the columns, and the genes by the rows. In the inverse case, the user can indicate this by giving the parameter `cellsOnCol = FALSE` to the function. *Note that in the Scala implementation, the samples (here cells) are represented by the rows, which is traditionally the case in machine learning.*


### Rareness score
To operate on large instances, MicroCellClust 2 uses a rareness score as heuristics. One should compute a vector `rs` that contains a value for each cell. Existing methods include :
* FiRE [Jindal *et al.*, 2018] https://github.com/princethewinner/FiRE
* DoRC [Chen *et al.*, 2019] https://github.com/chenxofhit/DoRC


## Running MicroCellClust

MicroCellClust 2 is run using the following function:
``` R
mcc.res = runMCC(mcc.rs, mat.filt, gene.sum = gs, rareness.score = rs) # If needed, don't forget to set the `data.tfo` and `cellesOnCol` arguments
```
As this function calls the Scala solver, the user must provide the `mcc.rs` object as first argument. `mcc.res` returns a list containing the indices and names of the cells (`mcc.res$cells.idx` and `mcc.res$cells.names`) and genes (`mcc.res$genes.idx` and `mcc.res$genes.names`) composing the identified bicluster, as well as the latter's objective value (`mcc.res$obj.value`). `mcc.res$info` contains additional information, including intermediate results.

### Optimization problem parameters

Arguments `kappa` and `nNeg` correspond to the $\kappa \ge 0$ and $\mu \in [0, 1]$ metaparameters of the optimisation problem [1], controlling respectively the out-of-cluster expression and the maximum proportion of negative values inside the bicluster. When omitted, default values of `kappa = 100 / ncol(mat)` and `nNeg = 0.1` are taken. 

These values are automatically tuned during the execution of the `runMCC` function. To disable this automatic tuning, set respectively `k.adapt = FALSE` and `n.adapt = FALSE`. To tune these values manually, one could re-run several times the `runMCC` function until obtaining the desired result. However, as it operates in two steps (beam search + local search) [2], one can perform this tuning by only re-running the second step (the local search) using:
``` R
runMCC.quick(mcc.rs, mat.filt, mcc.res$cells.idx, mcc.res$genes.idx, gene.sum = gs)
```

### Looking for other solutions

Once a first cell/gene bicluster has been found, one can re-run `runMCC` while excluding the cells in the `mcc.res` solution to search for another bicluster. This can be done using:
``` R
mcc.res.2 = runMCC(mcc.rs, mat.filt, gene.sum = gs, rareness.score = rs, cells.excl = mcc.res$cells.idx) # If needed, don't forget to set the `data.tfo` and `cellesOnCol` arguments
```
Alternatively, one could also/instead exclude the genes from the first solution using `genes.excl = mcc.res$genes.idx`.