# MicroCellClust

This repository contains an implementation of the MicroCellClust optimisation problem presented in Bioinformatics, btab239. It is designed to search for rare cell-clusters (typically less than 10% of the number of cells) that present a highly specific gene expression in scRNA-seq data. It has been implemented in the Scala programming language. This repository also contains an interface to run MicroCellClust in R.

Version 1.2 is also available on Zenodo http://doi.org/10.5281/zenodo.4580332

### References

A. Gerniers, O. Bricard and P. Dupont (2021). MicroCellClust: mining rare and highly specific subpopulations from single-cell expression data. *Bioinformatics*, 37(19), 3220-3227. https://doi.org/10.1093/bioinformatics/btab239



## Using MicroCellClust in R

### RScala

To run MicroCellClust in R, the user must install the `rscala` package (https://cran.r-project.org/web/packages/rscala/index.html), that can call the Scala implementation of MicroCellClust from within R. Once installed, he can check wether there is a version of Scala on his computer using:

``` R
rscala::scalaConfig()
```
This function tries to find Scala and Java on the user's computer and, if needed, downloads and installs Scala and Java in the user's `~/.rscala` directory.


### Initializing MicroCellClust

To run MicroCellClust, the user must first set his working directory (`setwd`) to the location of the `MicroCellClust`directory, and then load the functions contained in the `micro_cell_clust.R` script.
``` R
source("./src/main/R/micro_cell_clust.R")
```
Then, he must initialise an instance of RScala containing the MicroCellClust source code.
``` R
mcc.rs = initMccRscala()
```
This `mcc.rs` object needs to be passed as the first argument to any R function using Scala code.


### Preprocessing the data matrix

Let's assume the expression data is contained in a dataframe called `exprData`. This data is assumed to contain positive values when a gene is expressed in a cell, and negative values if it is not expressed (a possible scaling of the (normalized) count expression data to obtain this is $\log_{10}(x + 0.1)$).

As discussed in the paper, the first step is removing genes that are expressed in too many cells of the dataset. The function `preprocessGenes` filters out genes that are expressed in more than a certain proportion `thresh` of the cells (typically 25%).
``` R
exprData.filt = preprocessGenes(exprData, thresh = 0.25)
```
**Remark:** By default, the R implementation represents the cells by the columns, and the genes by the rows. In the opposite case, the user can indicate this by giving the parameter `cellsOnCol = FALSE` to the function. *Note that in the Scala implementation, the samples (here cells) are represented by the rows, which is usually the case in machine learning.*


### Running MicroCellClust

Running MicroCellClust is done using the following function:
``` R
mcc.res = runMCC(mcc.rs, exprData.filt, kappa = 1, nNeg = 0.1)
```
As this function calls the Scala solver, the user must provide the `mcc.rs` object as first argument. Arguments `kappa` and `nNeg` correspond to the $\kappa$ and $\mu$ metaparameters of the optimisation problem, controlling respectively the out-of-cluster expression and the maximum proportion of negative values inside the bicluster. We suggest initial values of $\kappa = \frac{100}{nb. cells}$ and $\mu = 0.1$.

`mcc.res` returns a list containing the indices and names of the cells (`mcc.res$cells.idx` and `mcc.res$cells.names`) and genes (`mcc.res$genes.idx` and `mcc.res$genes.names`) composing the identified bicluster, as well as its objective value (`mcc.res$obj.value`).


#### Parameter for the heuristic solver
The optimisation problem described in the paper is a NP-hard problem. To solve it efficiently, we use a heuristic solver that drives the search to the most promising zones of the search space. The solver follows a breadth-first search strategy, meaning it first evaluates all possible biclusters composed of 2 cells (level 2), then biclusters composed of 3 cells (level 3), etc. At each level, the heuristic will select only a fraction of the biclusters, the ones with above average objective values, to continue the search. Since the distribution of objective values roughly follow a power law, this corresponds to ignoring the long tail of the distribution (see supplementary data of the paper).

To do so, the `nHeuristic` parameter must be set to the rank in the distribution roughly corresponding to the beginning of the long tail:
``` R
pairs.obj = estimateHeuristic(mcc.rs, exprData.filt, kappa = 1)
plot(1:100, pairs.obj[1:100]) 
# Choose a value roughly corresponding to the beginning of the long tail
mcc.res = runMCC(mcc.rs, exprData.filt, kappa = 1, nNeg = 0.1, nHeuristic = 20)
```

#### Additional parameters (not used in the paper)

* `stopNoImprove`(int, default `25`). By default, the solver stops the search when no improvement to the current best solution is found during 25 successive search levels. This can be changed setting the `stopNoImprove` parameter to another value.
* `maxNbCells`(int, default `Inf`). Alternatively, a fixed upper limit on the number of cells can be defined by setting `maxNbCells`. The search is stopped once this level is reached.
* `minCoExpGenes`(int, default `0`). Adds an extra constraint specifying a minimum number of genes that must be expressed in all cells of the cluster to form a valid solution. By default, this constraint is not present.
* `verbose`(boolean, default `TRUE`). Printing information during execution can be disabled by setting the argument `verbose = FALSE`.

#### Top-$k$ search

Once a bicluster has been found, MicroCellClust may be run again to search for another possible solution among the remaining cells:
``` R
previous.cells = mcc.res$cells.idx
mcc.res.2 = runMCC.exclude(mcc.rs, exprData.filt, excl = previous.cells, kappa = 1, nNeg = 0.1)
```

## Using MicroCellClust in Scala

MicroCellClust can be executed directly in Scala. **Note that for the Scala implementation, the data must be represented with cells (samples) as rows genes as columns.** Assuming `exprData` is a 2D-array of doubles, and `cellNames` and `geneNames` are arrays of strings:
``` Scala
import Preprocessor.expressionFilter
import Solver.findCluster

val exprDataFilt, geneNamesFilt = expressionFilter(exprData, geneNames, thresh = 0.25)

val resCellsID, resGenesID, resObj = findCluster(exprDataFilt, kappa = 1, nNeg = 0.1)

val resCellsNames = resCellsID.map(cellNames(_))
val resGenesNames = resGenesID.map(geneNamesFilt(_))
```
