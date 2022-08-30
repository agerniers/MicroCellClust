import Objective._

import scala.collection.mutable.{PriorityQueue, Set => MutSet}

/*
 * Author: Alexander Gerniers (UCLouvain)
 */
object Solver {

    /**
      * Find the bicluster maximising the objective using a heuristic search
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param preMarkSum the sum of positive values for each marker (required if m contains only part of the samples)
      * @param rarenessScore a list of rareness scores for each sample
      * @param nNeg the maximum percentage of negative values allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @param kAdapt whether to adapt kappa if too many/few markers are selected
      * @param nHeurPair the number of pairs to form for each sample at level 2
      * @param nHeurKeep the number of top-solutions to consider for expansion at the next level
      * @param nhAdapt whether to decrease nHeurKeep during the search
      * @param maxNbSam the maximum number of samples in the bicluster
      * @param stopNoImprove stop search when no better solution is found after x levels
      * @param verbose enable/disable printing
      * @return - an assignment of samples
      *         - an assignment of markers
      *         - the corresponding objective value
      */
    def findCluster(m: Array[Array[Double]], preMarkSum: Array[Double] = Array(), rarenessScore: Array[Double] = Array(),
                    nNeg: Double = 0.1, kappa: Double = 1, kAdapt: Boolean = true, nHeurPair: Int = 100, nHeurKeep: Int = 100,
                    nhAdapt: Boolean = true, maxNbSam: Int = Int.MaxValue, stopNoImprove: Int = 25, verbose: Boolean = true): (List[Int], List[Int], Double, Double) = {
        val expr = buildExprMap(m)
        val markSum = if (preMarkSum.length == m(0).length) preMarkSum else getMarkSum(m)
        val nSam = m.length
        val hOrd = rarenessScore.zipWithIndex.sortBy(-_._1)

        var nH = nHeurKeep
        var kA = kappa

        var check = false

        if (verbose) {
            val l = "Level"
            val k = "(kappa)"
            val n = "(nb. heur.)"
            val p = "(nb. cand.)"
            val nb = "New best?"
            val g = "Nb. markers"
            val o = "Obj. value"
            val m = "(mean sim.)"
            val t = "Exe. time [s]"
            println(f"| $l%-6s | $k%-8s | $n%-12s | $p%-12s | $nb%-12s | $g%-12s | $o%-12s | $m%-12s | $t%-15s |")
        }

        // Evaluation of pairs of samples
        val nBestQueue = PriorityQueue[(Set[Int], List[Int], Double)]()(Ordering[Double].on(x => -x._3))

        while (!check) {
            if (verbose) {
                val l = 2
                print(f"| $l%-6s | $kA%-8s | $nHeurPair%-12s | $nSam%-12s |")
            }
            val t0 = System.currentTimeMillis

            if (hOrd.length == nSam) { // Use rareness score
                for (a <- (0 until nSam); b <- ((a + 1) until (nSam min (a + nHeurPair)))) {
                    val p = Set(hOrd(a)._2, hOrd(b)._2)
                    val (markers, obj) = getMarkers(m, p.toList, expr, markSum, nNeg, kA)
                    nBestQueue += ((p, markers, obj))
                    if (nBestQueue.size > nH) nBestQueue.dequeue
                }
            } else { // Generate all pairs
                for (p <- (0 until nSam).combinations(2)) {
                    val (markers, obj) = getMarkers(m, p.toList, expr, markSum, nNeg, kA)
                    nBestQueue += ((p.toSet, markers, obj))
                    if (nBestQueue.size > nH) nBestQueue.dequeue
                }
            }

            val t1 = System.currentTimeMillis

            val cBest = nBestQueue.toList.reverse(0)

            if (verbose) {
                val nb = "    ****"
                val g = cBest._2.length
                val o = (cBest._3 * 100).round / 100.0
                val m = ""
                val t = (t1 - t0).toDouble / 1000
                println(f" $nb%-12s | $g%-12s | $o%-12s | $m%-12s | $t%-15s |")
            }

            // Check if enough markers are selected (otherwise decrease kappa temporarily)
            if (kAdapt && cBest._2.length < 20) {
                var nTry = 5
                var kTest = kA / 2
                while (nTry > 0) {
                    val ng = getMarkers(m, cBest._1.toList, expr, markSum, nNeg, kTest)._1.length
                    if (ng < 20) {
                        nTry -= 1
                        kTest = kTest / 2
                    } else if (ng > 100) {
                        nTry -= 1
                        kTest = kTest * 1.5
                    } else {
                        nTry = 0
                    }
                }
                kA = kTest
            } else {
                check = true
            }
        }

        val nBest = nBestQueue.toList.reverse
        var best = nBest(0)
        var prevLvlNBest = nBest.map(x => x._1)

        var lvl = 3
        val maxLvl = maxNbSam min nSam
        var noImprove = 0
        var finished = false
        var cellFilt = false
        var samPool = (0 until nSam).toSet
        while (!finished && lvl <= maxLvl) {
            if (verbose) {
                val p = samPool.size
                print(f"| $lvl%-6s | $kA%-8s | $nH%-12s | $p%-12s |")
            }
            val t0 = System.currentTimeMillis
            val nBestQueue = PriorityQueue[(Set[Int], List[Int], Double)]()(Ordering[Double].on(x => -x._3))

            val hist = MutSet[Set[Int]]()

            for (p <- prevLvlNBest.indices) {
                val prev = prevLvlNBest(p)
                val prevExpSums = getExpSums(m, prev.toList, nNeg)
                for (i <- samPool if !prev.contains(i)) {
                    val samples = prev + i
                    val l = hist.size
                    hist += samples
                    if (hist.size - l >= 1) {
                        val (markers, obj) = getMarkersFromPrev(m, samples.size, i, markSum, prevExpSums, nNeg, kA)
                        nBestQueue += ((samples, markers, obj))
                        if (nBestQueue.size > nH) nBestQueue.dequeue
                    }
                }
            }
            val t1 = System.currentTimeMillis

            if (nBestQueue.size > 0) {
                val nBest = nBestQueue.toList.reverse
                val cBest = nBest(0)


                val msim = if (nH <= 1) 1
                else (nBest.drop(1).map(x => ((x._2 intersect cBest._2).size.toDouble / cBest._2.size)).sum / (nBest.length - 1))

                // If similarity of markers between solution is high: decrease nHeurKeep
                if (nhAdapt) {
                    if (msim >= 0.98 && lvl >= 200) nH = 1
                    else if (msim >= 0.98) nH = (nHeurKeep / 8) max 1
                    else if (msim >= 0.9) nH = nHeurKeep / 4
                    else if (msim >= 0.7) nH = nHeurKeep / 2
                    else nH = nHeurKeep
                }

                // If similarity of markers between solution is high: keep only samples that express these markers
                val t2 = System.currentTimeMillis()

                if (!cellFilt && msim >= 0.9) {
                    cellFilt = true
                    val markUnion = nBest.map(_._2).flatten.distinct
                    samPool = (0 until nSam).toSet.filter(i => getSupport(i, markUnion, expr) >= 0.3)
                } else if (cellFilt && lvl % 10 == 0) {
                    val markUnion = nBest.map(_._2).flatten.distinct
                    samPool = (0 until nSam).toSet.filter(i => getSupport(i, markUnion, expr) >= 0.3)
                }
                val t3 = System.currentTimeMillis()

                if (cBest._3 > best._3) {
                    best = (cBest._1, cBest._2, cBest._3)
                    noImprove = 0
                    if (verbose) {
                        val nb = "    ****"
                        val g = cBest._2.length
                        val o = (cBest._3 * 100).round / 100.0
                        val m = (msim * 10000).round / 100.0
                        val t = (t1 - t0 + t3 - t2).toDouble / 1000
                        println(f" $nb%-12s | $g%-12s | $o%-12s | $m%-12s | $t%-15s |")
                    }

                    if (kAdapt) {
                        if (cBest._2.length <= 200) {
                            if (kA < kappa) { // If kappa was temporarily decreased, check if we can go back to normal
                                val ng = getMarkers(m, cBest._1.toList, expr, markSum, nNeg, kappa)._1.length

                                if (ng >= 20) {
                                    kA = kappa
                                    best = (Set(), List(), 0.0)
                                } else {
                                    prevLvlNBest = nBest.map(x => x._1)
                                    lvl += 1
                                }
                            } else {
                                prevLvlNBest = nBest.map(x => x._1)
                                lvl += 1
                            }
                        } else { // If the solution contains too many markers, increase kappa
                            var nTry = 5
                            var kTest = kA * 2
                            while (nTry > 0) {
                                val ng = getMarkers(m, cBest._1.toList, expr, markSum, nNeg, kTest)._1.length
                                if (ng > 100) {
                                    nTry -= 1
                                    kTest = kTest * 2
                                } else if (ng < 50) {
                                    nTry -= 1
                                    kTest = kTest * 0.75
                                } else {
                                    nTry = 0
                                }
                            }
                            kA = kTest
                            best = (Set(), List(), 0.0)
                        }
                    } else {
                        prevLvlNBest = nBest.map(x => x._1)
                        lvl += 1
                    }
                } else {
                    noImprove += 1
                    if (verbose) {
                        val nb = ""
                        val g = cBest._2.length
                        val o = (cBest._3 * 100).round / 100.0
                        val m = (msim * 10000).round / 100.0
                        val t = (t1 - t0 + t3 - t2).toDouble / 1000
                        println(f" $nb%-12s | $g%-12s | $o%-12s | $m%-12s | $t%-15s |")
                    }
                    if (noImprove >= stopNoImprove) {
                        finished = true
                        println()
                    } else if (cBest._3 <= 0.0) {
                        finished = true
                        println()
                    }
                    prevLvlNBest = nBest.map(x =>  x._1)
                    lvl += 1
                }
            } else {
                finished = true
                println()
            }
        }

        return (best._1.toList.sorted, best._2, best._3, kA)
    }

}
