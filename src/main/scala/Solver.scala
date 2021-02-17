import scala.collection.mutable.{ListBuffer, PriorityQueue}

import Objective.{buildExprMap, getCoExprMarks, getExpSums, getMarkSum, getMarkers, getMarkersFromPrev}

/*
 * Author: Alexander Gerniers
 */
object Solver {

    /**
      * Find the bicluster maximising the objective using a heuristic search
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param nNeg the maximum percentage of -1 allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @param nHeuristic the initial number of top-solutions to consider for expansion at the next level
      * @param maxNbSam the maximum number of samples in the bicluster
      * @param stopNoImprove stop search when no better solution is found after x levels
      * @param minCoExpMark a minimum number of markers that must be expressed in every row of the bicluster
      * @param excl a list of samples that need to be excluded from the search
      * @param verbose enable/disable printing
      * @return - an assignment of samples
      *         - an assignment of markers
      *         - the corresponding objective value
      */
    def findCluster(m: Array[Array[Double]], nNeg: Double = 0.1, kappa: Double = 1, nHeuristic: Int = 20,
                    maxNbSam: Int = Int.MaxValue, stopNoImprove: Int = 25, minCoExpMark: Int = 0,
                    excl: List[Int] = List(), verbose: Boolean = true): (List[Int], List[Int], Double) = {
        val expr = buildExprMap(m)
        val markSum = getMarkSum(m)
        val nSam = m.length
        val nSamMinusExcl = nSam - excl.length

        if (verbose) {
            println("Current search level: 2")
            println("\t " + ((nSamMinusExcl * nSamMinusExcl - nSamMinusExcl) / 2) + " pairs to evaluate")
        }
        val t0 = System.currentTimeMillis

        // Evaluation of all the pairs of cells
        val nBestQueue = PriorityQueue[(List[Int], List[Int], Double)]()(Ordering[Double].on(x => -x._3))
        for (p <- (0 until nSam).combinations(2)) {
            if (excl.length == 0 || !excl.exists(p contains _)) {
                val (markers, obj) = getMarkers(m, p.toList, expr, markSum, nNeg, kappa)
                nBestQueue += ((p.toList, markers, obj))
                if (nBestQueue.size > nHeuristic) nBestQueue.dequeue
            }
        }

        val nBest = nBestQueue.dequeueAll.toList.reverse
        var best = nBest(0)
        var prevLvlNBest = nBest.map(x =>  (x._1, getCoExprMarks(x._1, expr)))
        val t1 = System.currentTimeMillis

        if (verbose) {
            println("\t NEW BEST: obj. val.: " + best._3)
            println("\t samples (idx. from 1): " + best._1.map(_ + 1).mkString(", "))
            println("\t nb. markers: " + best._2.length)
            println("\t Computation time [s]: " + ((t1 - t0).toDouble / 1000))
        }

        var lvl = 3
        val maxLvl = maxNbSam min nSam
        var noImprove = 0
        var finished = false
        while (!finished && lvl <= maxLvl) {
            if (verbose) {
                println("Current search level: " + lvl)
            }
            val t0 = System.currentTimeMillis
            val nBestQueue = PriorityQueue[(List[Int], List[Int], Double, List[Int])]()(Ordering[Double].on(x => -x._3))

            for (p <- prevLvlNBest.indices) {
                val prevExpSums = getExpSums(m, prevLvlNBest(p)._1, nNeg)
                var toExclude = excl ++ prevLvlNBest(p)._1
                for (pp <- 0 until p) {
                    val setDif = prevLvlNBest(pp)._1.filterNot(prevLvlNBest(p)._1.toSet)
                    if (setDif.length == 1) toExclude ++= setDif
                }
                for (i <- 0 until nSam if !toExclude.contains(i)) {
                    val samples = (i :: prevLvlNBest(p)._1).sorted
                    val coExp = (prevLvlNBest(p)._2.toSet intersect expr(i).toSet).toList
                    if (coExp.length >= minCoExpMark) {
                        val (markers, obj) = getMarkersFromPrev(m, samples, i, markSum, prevExpSums, nNeg, kappa)
                        nBestQueue += ((samples, markers, obj, coExp))
                        if (nBestQueue.size > nHeuristic) nBestQueue.dequeue
                    }
                }
            }

            if (nBestQueue.size > 0) {
                val nBest = nBestQueue.dequeueAll.toList.reverse
                prevLvlNBest = nBest.map(x => (x._1, x._4))
                if (nBest(0)._3 > best._3) {
                    best = (nBest(0)._1, nBest(0)._2, nBest(0)._3)
                    noImprove = 0
                    if (verbose) {
                        println("\t NEW BEST: obj. val.: " + best._3)
                        println("\t samples (idx. from 1): " + best._1.map(_ + 1).mkString(", "))
                        println("\t nb. markers: " + best._2.length)
                    }
                } else {
                    noImprove += 1
                    if (noImprove >= stopNoImprove) {
                        finished = true
                        if (verbose) {
                            println("\t No improvement after " + stopNoImprove + " levels: search stopped")
                        }
                    }
                }
                lvl += 1
            } else {
                finished = true
            }

            val t1 = System.currentTimeMillis
            if (verbose) {
                println("\t Computation time [s]: " + ((t1 - t0).toDouble / 1000))
            }
        }

        return best
    }

    /**
      * Get the objective value of each pair of samples
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return the objective values in decreasing order
      */
    def evaluatePairs(m: Array[Array[Double]], kappa: Double = 1): List[(List[Int], Double)] = {
        val expr = buildExprMap(m)
        val markSum = getMarkSum(m)
        val nSam = m.length

        val lb = ListBuffer[(List[Int], Double)]()

        for (p <- (0 until nSam).combinations(2)) {
            val (_, obj) = getMarkers(m, p.toList, expr, markSum, 0, kappa)
            lb += ((p.toList, obj))
        }

        return lb.toList.sortBy(-_._2)
    }

}
