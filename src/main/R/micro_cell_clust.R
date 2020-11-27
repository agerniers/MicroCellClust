# Author: Alexander Gerniers
library("rscala")


#' Initialises the RScala instance to run MicroCellClust.
#' 
#' @return An instance of RScala for MicroCellClust
#' 
initMccRscala = function() {
  mcc.rs = scala()
  
  processFile = function(filepath) {
    f = file(filepath, "r")
    lines = readLines(f, n = -1)
    paste(lines, collapse = "\n")
  }
  
  mcc.rs + processFile("./src/main/scala/Objective.scala")
  mcc.rs + processFile("./src/main/scala/Solver.scala")
  mcc.rs
}


#' Evaluates the objective values of each pair of cells and returns them in decreasing order. 
#' This can be used to estimate an initial number for the heuristic.
#' 
#' @param mcc.rs A MicroCellClust RScala instance
#' @param data An expression matrix
#' @param cellsOnCol Wether the cells are described by the columns (\code{TRUE}, default) or the rows (\code{FALSE})
#' @param kappa The out-of-cluster expression penalty constant
#' 
#' @return A vector of objective values of the pairs in decreasing order
#' 
estimateHeuristic = function(mcc.rs, data, cellsOnCol = TRUE, kappa = 1.0) {
  data.matrix = as.matrix(data)
  if (cellsOnCol) {
    data.matrix = t(data.matrix)
  }
  
  res = mcc.rs(m = data.matrix, k = kappa) * 'Solver.evaluatePairs(m, k)'
  
  pairs.id = mcc.rs(res = res) * 'res.map(_._1.mkString("-")).toArray'
  pairs.obj = mcc.rs(res = res) * 'res.map(_._2).toArray'
  names(pairs.obj) = pairs.id
  pairs.obj
}


#' Runs MicroCellClust on the given data to find the bicluster maximising the objective using heuristic search
#' 
#' @param mcc.rs A MicroCluster RScala instance
#' @param data An expression matrix
#' @param cellsOnCol Wether the cells are described by the columns (\code{TRUE}, default) or the rows (\code{FALSE})
#' @param kappa The out-of-cluster expression penalty constant (default: 1)
#' @param nNeg The maximum percentage of -1 allowed inside the cluster (default: 0.1)
#' @param maxNbCells The maximum number of cells allowed to form a solution (past this number, the seach is stopped; default: Inf)
#' @param minCoExpGenes Requires a solution to contain at least a certain number of genes expressed in every cell of the bicluster (default: 0)
#' @param nbHeurInit Parameter for the heuristic: the initial number of top-solutions to consider for expansion at the next level
#' @param nbHeurAdd Parameter for the heuristic: the number of extra top-solutions to consider at each level
#' 
#' @return The cells and genes forming the bicluster, and its objective value
#' 
runMCC = function(mcc.rs, data, cellsOnCol = TRUE, kappa = 1.0, nNeg = 0.1, maxNbCells = Inf, minCoExpGenes = 0, nbHeurInit = 20, nbHeurAdd = 10) {
  data.matrix = as.matrix(data)
  if (cellsOnCol) {  # Scala implementation considers cells on rows and genes on columns
    data.matrix = t(data.matrix)
  }
  
  res = mcc.rs(m = data.matrix, k = kappa, n = nNeg, mc = maxNbCells, mg = minCoExpGenes, nh = nbHeurInit, dh = nbHeurAdd) *
    'Solver.findCluster(m, kappa = k, nNeg = n, maxNbSam = mc.toInt, minCoExpMark = mg.toInt, nHeuristic = nh.toInt, deltaHeuristic = dh.toInt)'
  
  clust.cells = mcc.rs(res = res) * 'res._1.map(_ + 1).toArray'
  clust.genes = mcc.rs(res = res) * 'res._2.map(_ + 1).toArray'
  clust.obj = mcc.rs(res = res) * 'res._3'
  
  if (cellsOnCol) {
    res.lst = list(clust.cells, colnames(data)[clust.cells], clust.genes, rownames(data)[clust.genes], clust.obj)
    names(res.lst) = c("cells.idx", "cells.names", "genes.idx", "genes.names", "obj.value")
    res.lst
  } else {
    res.lst = list(clust.cells, rownames(data)[clust.cells], clust.genes, colnames(data)[clust.genes], clust.obj)
    names(res.lst) = c("cells.idx", "cells.names", "genes.idx", "genes.names", "obj.value")
    res.lst
  }
}

#' Get the optimal assignment of genes for a given cell-cluster
#' 
#' @param mcc.rs A MicroCluster RScala instance
#' @param data An expression matrix
#' @param cells.idx A list of cells
#' @param cellsOnCol Wether the cells are described by the columns (\code{TRUE}, default) or the rows (\code{FALSE})
#' @param kappa The out-of-cluster expression penalty constant (default: 1)
#' @param nNeg The maximum percentage of -1 allowed inside the cluster (default: 0.1)
#' 
#' @return The correspoonding cluster of genes and the objective value
#' 
getGenes = function(mcc.rs, data, cells.idx, cellsOnCol = TRUE, kappa = 1.0, nNeg = 0.1) {
  data.matrix = as.matrix(data)
  if (cellsOnCol) {
    data.matrix = t(data.matrix)
  }
  
  res = mcc.rs(m = data.matrix, c = cells.idx, k = kappa, n = nNeg) *
    'Objective.getMarkers(m, c.map(_.toInt - 1).toList, Objective.buildExprMap(m), Objective.getMarkSum(m), kappa = k, nNeg = n)'
  
  clust.genes = mcc.rs(res = res) * 'res._1.map(_ + 1).toArray'
  clust.obj = mcc.rs(res = res) * 'res._2'
  
  if (cellsOnCol) {
    res.lst = list(clust.genes, rownames(data)[clust.genes], clust.obj)
    names(res.lst) = c("genes.idx", "genes.names", "obj.value")
    res.lst
  } else {
    res.lst = list(clust.genes, colnames(data)[clust.genes], clust.obj)
    names(res.lst) = c("genes.idx", "genes.names", "obj.value")
    res.lst
  }
}

#' Preprocesses the data matrix to remove genes that are expressed in too many cells. 
#' Filters out the genes that are expressed in only 0 or 1 cell, 
#' and those expressed in more than a certain proportion \code{thresh} of the cells.
#' 
#' @param data An expression matrix
#' @param cellsOnCol Wether the cells are described by the columns (\code{TRUE}, default) of the rows (\code{FALSE})
#' @param type The type of filtering that is used
#' @param thresh The threshold value: maximal proportion of cells in which a gene may be expressed
#' 
#' @return The preprocessed data matrix
#' 
preprocessGenes = function(data, cellsOnCol = TRUE, type = "remove", thresh = 0.25) {
  nb_expr = function(x) {
    length(which(x >= 0))
  }
  if (type == "remove") {
    if (cellsOnCol) {
      data.expr = apply(data, 1, nb_expr)
      data[which(data.expr >= 2 & data.expr <= thresh * ncol(data)), ]
    } else {
      data.expr = apply(data, 2, nb_expr)
      data[, which(data.expr >= 2 & data.expr <= thresh * nrow(data))]
    }
  } else {
    data
  }
}
