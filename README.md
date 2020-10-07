# MicroCellClust

This repository contains an implementation of the MicroCellClust optimisation problem presented in [...]. It was designed to search for small cell-clusters (typically less than 10% of the number of cells) that present a highly specific gene expression in scRNA-seq data. It has been implemented in the Scala programming language. This repository also contains an interface to run MicroCellClust in R.



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

Let's assume the expression data is contained in a dataframe called `exprData`. As discussed in the paper, the first step is to get rid of genes that are expressed in too many cells of the dataset. The function `preprocessGenes` filters out genes that are expressed in more than a certain proportion `thresh` of the cells.
``` R
exprData.filt = preprocessGenes(exprData, thresh = 0.25)
```
**Remark:** By default, the R implementation represents the cells by the columns, and the genes by the rows. In the inverse case, the user can indicate this by giving the parameter `cellsOnCol = FALSE` to the function. *Note that in the Scala implementation, the samples (here cells) are represented by the rows, which is traditionally the case in machine learning.*


### Running MicroCellClust

MicroCellClust is run using the following function:
``` R
mcc.res = runMCC(mcc.rs, exprData.filt, kappa = 1, nNeg = 0.1)
```
As this function calls the Scala solver, the user must provide the `mcc.rs` object as first argument. Arguments `kappa` and `nNeg` correspond to the $\kappa$ and $\mu$ methaparameters of the optimisation problem, controlling respectively the out-of-cluster expression and the maximum proportion of negative values inside the bicluster. `mcc.res` returns a list containing the indices and names of the cells (`mcc.res$cells.idx` and `mcc.res$cells.names`) and genes (`mcc.res$genes.idx` and `mcc.res$genes.names`) composing the identified bicluster, as well as the latter's objective value (`mcc.res$obj.value`).

#### Additional parameters
Additional parameters can be used. A maximum number of cells that the bicluster may contain can be enforced using the `maxNbCells` argument. Once the search reaches clusters of this size, it will stop. Since the goal of MicroCellClust is to find relatively small clusters, this parameter is useful to avoid exploring parts of the search space that are anyway uninteristing.
``` R
mcc.res = runMCC(mcc.rs, exprData.filt, kappa = 1, nNeg = 0.1, maxNbCells = 50)
```
An other possibly useful parameter is `minCoExpGenes`. It sets a minimum number of genes that must be expressed in every cell of the bicluster to form a solution. Since we search for highly specific clusters, it can make sense to impose that a bicluster contains at least a few genes without any negative expression. This constraint has the advantage that it will eliminate a lot of biclusters during the search, which can give a significant speedup.
``` R
mcc.res = runMCC(mcc.rs, exprData.filt, kappa = 1, nNeg = 0.1, minCoExpGenes = 2)
```

### Parameters for the Heuristic
The optimisation problem described in the paper is a NP-hard problem. To solve it efficiently, we use a heuristic that drives the search to the most promissing zones of the search space. The solver follows a breadth-first search strategy, meaning it first evaluates all possible biclusters composed of 2 cells (level 2), then biclusters composed of 3 cells (level 3), etc. At each level, the heuristic will select only a fraction of the biclusters, the ones with the highest objective values, to continue the search. Initially, a number `nbHeurInit` of solutions are selected, and at each subsequent level, this number is increased by `nbHeurAdd`:

* At level 2, all the pairs of cells are evaluated. Then, the `nbHeurInit` pairs with the highest objective value are selected.
* At level 3, biclusters are created by adding 1 cell to the biclusters that where selected at the previous level. These biclusters are evaluated and the `nbHeurInit` + `nbHeurAdd` best ones are selected.
* At level 4, we again add 1 cell to the previously selected biclusters, and select the `nbHeurInit`+ 2 * `nbHeurAdd` best ones among them.
* Etc.

The values `nbHeurInit` and `nbHeurAdd` will influence the duration of the search. The smaller they are, the faster the search will end as fewer possible solutions are evaluated. Of course, these values must not be too small, otherwise we might miss the optimal solution. A tradeoff should be found between searching in a sufficiently large part of the search space, and ensuring the search doesn't take too much time.

In order to have an estimation of `nbHeurInit`, the user can plot the objective values of all the pairs of cells, ordered decreasingly, and look at the evolution of this curve. This curve generally roughly follows a power law. A good estimation is to take the index where the slope becomes small as value for `nbHeurInit`.
``` R
pairs.obj = estimateHeuristic(mcc.rs, exprData.filt, kappa = 1)
plot(1:100, pairs.obj[1:100]) 
# For nbHeurInit, choose a value roughly corresponding to the place where the slope gets small
# For nbHeurAdd, choose for example nbInit/2
mcc.res = runMCC(mcc.rs, exprData.filt, kappa = 1, nNeg = 0.1, nbHeurInit = 20, nbHeurAdd = 10)
```