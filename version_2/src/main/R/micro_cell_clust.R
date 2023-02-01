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
    mcc.rs + processFile("./src/main/scala/LocalSearch.scala")
    mcc.rs
}




#' Computes the sum of positive expression for each gene in the transformed data
#'
#' @param data The (normalized) count data (will be transformed inside the function according to \code{data.tfo})
#' @param data.tfo How to transform the data in a matrix with positive and negative values. Either:
#'                 - "log10", i.e. log10(x + 0.1) (default)
#'                 - "log2", i.e. log2(x + 0.5)
#'                 - a function with one parameter (the original data) that returns the transformed data
#' @param cellsOnCol \code{TRUE} if data is GENES x CELLS (default), \code{FALSE} if CELLS x GENES
#'
#' @return A vector with values for each gene
#'
geneSum = function(data, data.tfo = "log10", cellsOnCol = TRUE) {
    if (data.tfo == "log10") {
        data.tfo = function (x) { log10(x + 0.1) }
    } else if (data.tfo == "log2") {
        data.tfo = function (x) { log2(x + 0.5) }
    } else if (data.tfo == "none") {
        data.tfo = function (x) { x }
    }
    count_cpe = function (x) {
        lx = data.tfo(x)
        sum(lx[lx > 0])
    }

    if (cellsOnCol) {
        if (ncol(data) < 50000) {
            as.vector(apply(data, 1, count_cpe))
        } else {
            r1 = as.vector(apply(data[, 1 : floor(ncol(data) / 2)], 1, count_cpe))
            r2 = as.vector(apply(data[, (floor(ncol(data) / 2) + 1) : ncol(data)], 1, count_cpe))
            r1 + r2
        }
    } else {
        if (nrow(data) < 50000) {
            as.vector(apply(data, 2, count_cpe))
        } else {
            r1 = as.vector(apply(data[1 : floor(nrow(data) / 2), ], 2, count_cpe))
            r2 = as.vector(apply(data[(floor(nrow(data) / 2) + 1) : nrow(data), ], 2, count_cpe))
            r1 + r2
        }
    }
}


#' Function calling the Scala code for the breadth-first solver
#'
#' @param mcc.rs A MicroCellClust RScala instance
#' @param data.matrix A CELLS x GENES expression matrix (with positive and negative values)
#' @param gene.sum The sum of positive expression in \code{data.matrix} for each gene (computed if not given)
#' @param rareness.score A list of rareness scores for each cell to use as search heuristic (optional)
#' @param nNeg The maximum proportion of negative values allowed inside the cluster (default: 0.1)
#' @param kappa The out-of-cluster expression penalty constant (default: 1)
#' @param k.adapt Whether to adapt kappa if too many/few markers are selected (default: \code{TRUE})
#' @param n.heur.pair The number of pairs to form for each sample at level 2 (default: 100)
#' @param n.heur.keep The number of top-solutions to consider for expansion at the next level (default: 100)
#' @param nh.adapt Whether to decrease \code{n.heur.keep} during the search (default: \code{TRUE})
#' @param max.nb.cells The maximum number of cells allowed to form a solution (past this number, the search is stopped) (default: Inf)
#' @param stop.no.improve stop search when no better solution is found after x levels (default: 25)
#' @param verbose enable/disable printing (default: \code{TRUE})
#'
#' @return The cells and genes forming the bicluster, its objective value, and the final value of kappa
#'
runMCC.bfs = function(mcc.rs, data.matrix, gene.sum, rareness.score = c(), nNeg = 0.1, kappa = 1.0, k.adapt = TRUE, n.heur.pair = 100, n.heur.keep = 100, nh.adapt = TRUE, max.nb.cells = Inf, stop.no.improve = 25, verbose = TRUE) {

    if (length(rareness.score) == nrow(data.matrix)) {
        res = mcc.rs(m = data.matrix, gs = gene.sum, rs = rareness.score, n = nNeg, k = kappa, ka = k.adapt, nhp = n.heur.pair, nhk = n.heur.keep, nha = nh.adapt, mc = max.nb.cells, sni = stop.no.improve, v = verbose) *
        'Solver.findCluster(m, preMarkSum = gs.toArray, rarenessScore = rs.toArray, nNeg = n, kappa = k, kAdapt = ka, nHeurPair = nhp.toInt, nHeurKeep = nhk.toInt, nhAdapt = nha, maxNbSam = mc.toInt, stopNoImprove = sni.toInt, verbose = v)'
    } else {
        res = mcc.rs(m = data.matrix, gs = gene.sum, n = nNeg, k = kappa, ka = k.adapt, nhp = n.heur.pair, nhk = n.heur.keep, nha = nh.adapt, mc = max.nb.cells, sni = stop.no.improve, v = verbose) *
        'Solver.findCluster(m, preMarkSum = gs.toArray, nNeg = n, kappa = k, kAdapt = ka, nHeurPair = nhp.toInt, nHeurKeep = nhk.toInt, nhAdapt = nha, maxNbSam = mc.toInt, stopNoImprove = sni.toInt, verbose = v)'
    }

    clust.cells = mcc.rs(res = res) * 'res._1.map(_ + 1).toArray'
    clust.genes = mcc.rs(res = res) * 'res._2.map(_ + 1).toArray'
    clust.obj = mcc.rs(res = res) * 'res._3'
    kappa.final = mcc.rs(res = res) * 'res._4'

    res.lst = list(clust.cells, rownames(data.matrix)[clust.cells], clust.genes, colnames(data.matrix)[clust.genes], clust.obj, kappa.final)
    names(res.lst) = c("cells.idx", "cells.names", "genes.idx", "genes.names", "obj.value", "kappa.final")
    res.lst
}


#' Function calling the Scala code for the local-search solver
#'
#' @param mcc.rs A MicroCellClust RScala instance
#' @param data.matrix A CELLS x GENES expression matrix (with positive and negative values)
#' @param init.sol A solution from which to start the local search (typically the result of \code{runMCC.bfs})
#' @param gene.sum The sum of positive expression in \code{data.matrix} for each gene (computed if not given)
#' @param nNeg The maximum proportion of negative values allowed inside the cluster (default: 0.1)
#' @param kappa The out-of-cluster expression penalty constant (default: 1)
#' @param nb.iter The number of iterations of each restart of the local search (default: 1000)
#' @param nb.restart The number of times the local search should restart from the initial solution (default: 10)
#' @param seed A seed for random events (default: 0, i.e. no seed is set)
#' @param verbose enable/disable printing (default: \code{TRUE})
#'
#' @return The cells and genes forming the bicluster, and its objective value
#'
runMCC.local = function(mcc.rs, data.matrix, init.sol, gene.sum, kappa = 1.0, nNeg = 0.1, nb.iter = 1000, nb.restart = 20, seed = 0, verbose = TRUE) {

    res = mcc.rs(m = data.matrix, is = init.sol, gs = gene.sum,  k = kappa, n = nNeg, ni = nb.iter, nr = nb.restart, s = seed, v = verbose) *
    'LocalSearch.localSearch(m, initSamples = is.map(_ - 1).toSet, preMarkSum = gs.toArray, kappa = k, nNeg = n, nbIter = ni.toInt, nbRestart = nr.toInt, seed = s.toInt, verbose = v)'

    clust.cells = mcc.rs(res = res) * 'res._1.map(_ + 1).toArray'
    clust.genes = mcc.rs(res = res) * 'res._2.map(_ + 1).toArray'
    clust.obj = mcc.rs(res = res) * 'res._3'

    res.lst = list(clust.cells, rownames(data.matrix)[clust.cells], clust.genes, colnames(data.matrix)[clust.genes], clust.obj)
    names(res.lst) = c("cells.idx", "cells.names", "genes.idx", "genes.names", "obj.value")
    res.lst
}


#' Run MicroCellClust
#'
#' @param mcc.rs A MicroCellClust RScala instance
#' @param data The (normalized) count data (i.e. positive values; will be transformed inside the function according to \code{data.tfo})
#'             If data already contains positive and negative values, set data.tfo to "none"
#' @param gene.sum The sum of positive expression in the transformed data for each gene (computed if not given)
#' @param rareness.score A list of rareness scores for each cell to use as search heuristic for
#'                       the breadth-first search (optional, but recommended for data with > 1000 cells)
#' @param data.tfo How to transform the data in a matrix with positive and negative values. Either:
#'                 - "log10", i.e. log10(x + 0.1) (default)
#'                 - "log2", i.e. log2(x + 0.5)
#'                 - a function with one parameter (the original data) that returns the transformed data
#'                 - "none" if the data has already been transformed
#' @param cells.excl Cells to exclude (e.g. if they are part of a previous solution)
#' @param genes.excl Genes to exclude (e.g. if they are part of a previous solution)
#' @param nNeg The maximum proportion of negative values allowed inside the cluster (default: 0.1)
#' @param k.adapt Whether to adapt nNeg before local search (default: \code{TRUE})
#' @param kappa The out-of-cluster expression penalty constant (default: "auto", i.e. computed as 100 / nb. cells)
#' @param k.adapt Whether to adapt kappa if too many/few markers are selected during breadth-first search (default: \code{TRUE})
#' @param iqr If nb.cells > 1000 and \code{rareness.score} is given, define the threshold of cell selection based on the rareness
#'            score (for breadth-first search) as Q3 + \code{iqr} * IQR (default: 1.5)
#' @param seed A seed for random events during the local search (default: 0, i.e. no seed is set)
#' @param verbose enable/disable printing (default: \code{TRUE})
runMCC = function(mcc.rs, data, gene.sum = c(), rareness.score = c(), data.tfo = "log10", cells.excl = c(), genes.excl = c(), nNeg = 0.1, n.adapt = TRUE, kappa = "auto", k.adapt = TRUE, iqr = 1.5, cellsOnCol = TRUE, seed = 0, verbose = TRUE) {
    if (kappa == "auto"  && cellsOnCol) {
        kappa = round(100 / ncol(data), 3)
    } else if (kappa == "auto" && !cellsOnCol) {
        kappa = round(100 / nrow(data), 3)
    }

    if (length(gene.sum) != nrow(data) && length(gene.sum) != ncol(data)) {
        gene.sum = geneSum(data, data.tfo, cellsOnCol)
    }
    if (data.tfo == "log10") {
        data.tfo = function (x) { log10(x + 0.1) }
    } else if (data.tfo == "log2") {
        data.tfo = function (x) { log2(x + 0.5) }
    } else if (data.tfo == "none") {
        data.tfo = function (x) { x }
    }

    # cell selection using rareness score
    if (length(rareness.score) > 1000) {
        th = quantile(rareness.score, 0.75) + iqr * IQR(rareness.score)
        cells.idx.1st = as.vector(setdiff(which(rareness.score >= th), cells.excl))
    } else if (cellsOnCol) {
        cells.idx.1st = as.vector(setdiff(1 : ncol(data), cells.excl))
    } else {
        cells.idx.1st = as.vector(setdiff(1 : nrow(data), cells.excl))
    }

    # remove genes not expressed in these cells
    if (data.tfo(42) == 42) { # I.e. data.tfo was "none"
        if (cellsOnCol) {
            genes.expr = rowSums(data[, cells.idx.1st] >= 0)
        } else {
            genes.expr = colSums(data[cells.idx.1st, ] >= 0)
        }
        genes.idx.1st = as.vector(setdiff(which(genes.expr > 2), genes.excl))
    } else { # Data assumed to contain only positive values
        if (cellsOnCol) {
            genes.expr = rowSums(data[, cells.idx.1st] > 0)
        } else {
            genes.expr = colSums(data[cells.idx.1st, ] > 0)
        }
        genes.idx.1st = as.vector(setdiff(which(genes.expr > 2), genes.excl))
    }

    if (verbose) {
        message(sprintf("Beginning 1st run with %d cells and %d genes", length(cells.idx.1st), length(genes.idx.1st)))
    }

    # 1st search: breadth-first
    if (cellsOnCol) {
        res.1st = runMCC.bfs(mcc.rs, data.tfo(as.matrix(t(data[genes.idx.1st, cells.idx.1st]))), gene.sum[genes.idx.1st], rareness.score[cells.idx.1st], nNeg = nNeg, kappa = kappa, k.adapt = k.adapt, verbose = verbose)
    } else {
        res.1st = runMCC.bfs(mcc.rs, data.tfo(as.matrix(data[cells.idx.1st, genes.idx.1st])), gene.sum[genes.idx.1st], rareness.score[cells.idx.1st], nNeg = nNeg, kappa = kappa, k.adapt = k.adapt, verbose = verbose)
    }

    # selection of interesting cells for 2nd search
    if (cellsOnCol) {
        cells.expr = colSums(data[genes.idx.1st[res.1st$genes.idx], ] > 0)
    } else {
        cells.expr = rowSums(data[, genes.idx.1st[res.1st$genes.idx]] > 0)
    }
    expr.thresh = max(0.4 * length(res.1st$genes.idx), min(0.5 * length(res.1st$genes.idx), min(cells.expr[cells.idx.1st[res.1st$cells.idx]])))
    # expr.thresh = 0.8 * min(cells.expr[cells.idx.1st[res.1st$cells.idx]])
    if (verbose) {
        message(sprintf("Selecting cells expressing %d percent of provided genes", round((expr.thresh / length(res.1st$genes.idx)) * 100)))
    }
    cells.idx.2nd = union(as.vector(setdiff(which(cells.expr >= expr.thresh), cells.excl)), as.vector(cells.idx.1st[res.1st$cells.idx]))

    # tuning of nNeg parameter
    if (n.adapt) {
        #neg.cand = intersect(as.vector(which(cells.expr >= 0.75 * length(res.1st$genes.idx))), cells.idx.1st[-res.1st$cells.idx])
        neg.cand = setdiff(as.vector(which(cells.expr >= 0.75 * length(res.1st$genes.idx))), cells.excl)
        if (cellsOnCol) {
            neg.mat = data.tfo(data[genes.idx.1st[res.1st$genes.idx], c(cells.idx.1st[res.1st$cells.idx], neg.cand)])
            neg.pc = apply(neg.mat, 1, function (x) {
                length(which(x < 0)) / length(x)
            })
        } else {
            neg.mat = data.tfo(data[c(cells.idx.1st[res.1st$cells.idx], neg.cand), genes.idx.1st[res.1st$genes.idx]])
            neg.pc = apply(neg.mat, 2, function (x) {
                length(which(x < 0)) / length(x)
            })
        }
        new.nNeg = max(round(neg.pc * 20) / 20)
        if (verbose & new.nNeg != nNeg) {
            message(sprintf("Changing value of nNeg to %f", new.nNeg))
        }
    } else {
        new.nNeg = nNeg
    }

    # remove genes not expressed in a majority of cells of previous solution
    if (cellsOnCol) {
        genes.expr = rowSums(data[, cells.idx.2nd] > 0)
    } else {
        genes.expr = colSums(data[cells.idx.2nd, ] > 0)
    }
    genes.idx.2nd = as.vector(setdiff(which(genes.expr > (0.9 - new.nNeg) * length(res.1st$cells.idx)), genes.excl))

    if (verbose) {
        message(sprintf("Beginning 2nd run with %d cells and %d genes", length(cells.idx.2nd), length(genes.idx.2nd)))
    }

    # 2nd search: local search
    prev.sol = which(cells.idx.2nd %in% as.vector(cells.idx.1st[res.1st$cells.idx]))
    if (cellsOnCol) {
        res.2nd = runMCC.local(mcc.rs, data.tfo(as.matrix(t(data[genes.idx.2nd, cells.idx.2nd]))), prev.sol, gene.sum[genes.idx.2nd], kappa = res.1st$kappa.final, nNeg = new.nNeg, nb.iter = length(cells.idx.2nd) * 100, nb.restart = 20, seed = seed, verbose = verbose)
    } else {
        res.2nd = runMCC.local(mcc.rs, data.tfo(as.matrix(data[cells.idx.2nd, genes.idx.2nd])), prev.sol, gene.sum[genes.idx.2nd], kappa = res.1st$kappa.final, nNeg = new.nNeg, nb.iter = length(cells.idx.2nd) * 100, nb.restart = 20, seed = seed, verbose = verbose)
    }
    info.lst = list(cells.idx.1st, genes.idx.1st, cells.idx.1st[res.1st$cells.idx], res.1st$cells.names, genes.idx.1st[res.1st$genes.idx], res.1st$genes.names, res.1st$obj.val, cells.idx.2nd, genes.idx.2nd, res.1st$kappa.final, new.nNeg)
    names(info.lst) = c("cells.idx.1st", "genes.idx.1st", "res1.cells.idx", "res1.cells.names", "res1.genes.idx", "res1.genes.names", "res1.obj.value", "cells.idx.2nd", "genes.idx.2nd", "kappa.final", "nNeg.final")

    res.lst = list(cells.idx.2nd[res.2nd$cells.idx], res.2nd$cells.names, genes.idx.2nd[res.2nd$genes.idx], res.2nd$genes.names, res.2nd$obj.val, info.lst)
    names(res.lst) = c("cells.idx", "cells.names", "genes.idx", "genes.names", "obj.value", "info")
    res.lst
}


runMCC.quick = function(mcc.rs, data, init.cells, init.genes, gene.sum = c(), data.tfo = "log10", cells.excl = c(), genes.excl = c(), nNeg = 0.1, kappa = "auto", cellsOnCol = TRUE, nb.restart = 20, seed = 0, verbose = TRUE) {
    if (kappa == "auto"  && cellsOnCol) {
        kappa = round(100 / ncol(data), 3)
    } else if (kappa == "auto" && !cellsOnCol) {
        kappa = round(100 / nrow(data), 3)
    }

    if (length(gene.sum) != nrow(data) && length(gene.sum) != ncol(data)) {
        gene.sum = geneSum(data, data.tfo, cellsOnCol)
    }
    if (data.tfo == "log10") {
        data.tfo = function (x) { log10(x + 0.1) }
    } else if (data.tfo == "log2") {
        data.tfo = function (x) { log2(x + 0.5) }
    } else if (data.tfo == "none") {
        data.tfo = function (x) { x }
    }

    # selection of interesting cells
    if (cellsOnCol) {
        cells.expr = colSums(data[init.genes, ] > 0)
    } else {
        cells.expr = rowSums(data[, init.genes] > 0)
    }
    expr.thresh = max(0.4 * length(init.genes), min(0.5 * length(init.genes), min(cells.expr[init.cells])))
    if (verbose) {
        message(sprintf("Selecting cells expressing %d percent of genes returned by 1st run", round((expr.thresh / length(init.genes)) * 100)))
    }
    cells.idx = union(as.vector(setdiff(which(cells.expr >= expr.thresh), cells.excl)), as.vector(init.cells))

    # remove genes not expressed in a majority of cells of previous solution
    if (cellsOnCol) {
        genes.expr = rowSums(data[, cells.idx] > 0)
    } else {
        genes.expr = colSums(data[cells.idx, ] > 0)
    }
    genes.idx = as.vector(setdiff(which(genes.expr > (0.9 - nNeg) * length(init.cells)), genes.excl))

    if (verbose) {
        message(sprintf("Beginning 2nd run with %d cells and %d genes", length(cells.idx), length(genes.idx)))
    }

    # local search
    prev.sol = which(cells.idx %in% as.vector(init.cells))
    if (cellsOnCol) {
        res = runMCC.local(mcc.rs, data.tfo(as.matrix(t(data[genes.idx, cells.idx]))), prev.sol, gene.sum[genes.idx], kappa = kappa, nNeg = nNeg, nb.iter = length(cells.idx) * 100, nb.restart = nb.restart, seed = seed, verbose = verbose)
    } else {
        res = runMCC.local(mcc.rs, data.tfo(as.matrix(data[cells.idx, genes.idx])), prev.sol, gene.sum[genes.idx], kappa = kappa, nNeg = nNeg, nb.iter = length(cells.idx) * 100, nb.restart = nb.restart, seed = seed, verbose = verbose)
    }
    res.lst = list(cells.idx[res$cells.idx], res$cells.names, genes.idx[res$genes.idx], res$genes.names, res$obj.val)
    names(res.lst) = c("cells.idx", "cells.names", "genes.idx", "genes.names", "obj.value")
    res.lst
}