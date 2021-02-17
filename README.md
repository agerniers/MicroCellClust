# MicroCellClust

This repository contains an implementation of the MicroCellClust optimisation problem presented in [...]. It is designed to search for small cell-clusters (typically less than 10% of the number of cells) that present a highly specific gene expression in scRNA-seq data. It has been implemented in the Scala programming language. This repository also contains an interface to run MicroCellClust in R.



## Using MicroCellClust in R

### RScala

To run MicroCellClust in R, the user must install the `rscala` package (https://cran.r-project.org/web/packages/rscala/index.html), that can call the Scala implementation of MicroCellClust from within R. Once installed, he can check wether there is a version of Scala on his computer using:

``` R
rscala::scalaConfig()
```
This function tries to find Scala and Java on the user's computer and, if needed, downloads and installs Scala and Java in the user's `~/.rscala` directory.


### Initialising MicroCellClust

To run MicroCellClust, the user must the user must first set his working directory (`setwd`) to the location of the `MicroCellClust`directory, and then load the functions contained in the `micro_cell_clust.R` script.
``` R
source("./src/main/R/micro_cell_clust.R")
```
Then, he must initialise an instance of RScala containing the MicroCellClust source code.
``` R
mcc.rs = initMccRscala()
```
This `mcc.rs` object needs to be passed as the first argument to any R function using Scala code.


### Preprocessing the data matrix

Let's assume the expression data is contained in a dataframe called `exprData`. As discussed in the paper, the first step is removing genes that are expressed in too many cells of the dataset. The function `preprocessGenes` filters out genes that are expressed in more than a certain proportion `thresh` of the cells.
``` R
exprData.filt = preprocessGenes(exprData, thresh = 0.25)
```
**Remark:** By default, the R implementation represents the cells by the columns, and the genes by the rows. In the opposite case, the user can indicate this by giving the parameter `cellsOnCol = FALSE` to the function. *Note that in the Scala implementation, the samples (here cells) are represented by the rows, which is usually the case in machine learning.*


### Running MicroCellClust

MicroCellClust is run using the following function:
``` R
mcc.res = runMCC(mcc.rs, exprData.filt, kappa = 1, nNeg = 0.1)
```
As this function calls the Scala solver, the user must provide the `mcc.rs` object as first argument. Arguments `kappa` and `nNeg` correspond to the $\kappa$ and $\mu$ methaparameters of the optimisation problem, controlling respectively the out-of-cluster expression and the maximum proportion of negative values inside the bicluster. Good initial values are $\kappa = \frac{100}{nb. cells}$ and $\mu = 0.1$.

`mcc.res` returns a list containing the indices and names of the cells (`mcc.res$cells.idx` and `mcc.res$cells.names`) and genes (`mcc.res$genes.idx` and `mcc.res$genes.names`) composing the identified bicluster, as well as the latter's objective value (`mcc.res$obj.value`).

By default, the solver stops the search when no improvement to the current best solution is found during 25 successive search levels. This can be changed using the `stopNoImprove` parameter:
``` R
mcc.res = runMCC(mcc.rs, exprData.filt, kappa = 1, nNeg = 0.1, stopNoImprove = 25)
```
Printing information during the execution of `runMCC` can be disabled by setting the argument `verbose = FALSE`.


#### Parameter for the Heuristic
The optimisation problem described in the paper is a NP-hard problem. To solve it efficiently, we use a heuristic that drives the search to the most promissing zones of the search space. The solver follows a breadth-first search strategy, meaning it first evaluates all possible biclusters composed of 2 cells (level 2), then biclusters composed of 3 cells (level 3), etc. At each level, the heuristic will select only a fraction of the biclusters, the ones with the highest objective values, to continue the search. 

The number of solutions kept is defined by the `nHeuristic` parameter. Its value will influence the duration of the search after the evaluation of all pairs. However, it should be sufficiently big so as to avoid missing the optimum. In order to have an estimation of `nHeuristic`, the user can plot the objective values of all the pairs of cells, ordered decreasingly, and look at the evolution of this curve, which roughly follows a power law. A good estimation is to take the rank where the slope becomes small as value for `nHeuristic`.
``` R
pairs.obj = estimateHeuristic(mcc.rs, exprData.filt, kappa = 1)
plot(1:100, pairs.obj[1:100]) 
# Choose a value roughly corresponding to the place where the slope gets small
mcc.res = runMCC(mcc.rs, exprData.filt, kappa = 1, nNeg = 0.1, nHeuristic = 20)
```

#### Top-$k$ search

Once a bicluster has been found, MicroCellClust may be re-run to search for another possible solution among the remaining cells:
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