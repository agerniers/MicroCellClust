import Objective._

import scala.math.{ceil, exp, floor, pow}
import scala.util.Random

/*
 * Author: Alexander Gerniers (UCLouvain)
 */
object LocalSearch {

    /**
      * Local search to improve initial solution
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param initSamples an initial assignment of samples
      * @param nNeg the maximum percentage of negative values allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @param nbIter the number of neighbours to visit
      * @param preMarkSum the sum of positive values for each marker (required if m contains only part of the samples)
      * @param seed the seed of the pseudo-random number generator (0 == no seed)
      * @param verbose enable/disable printing
      * @return
      */
    def localSearch(m: Array[Array[Double]], initSamples: Set[Int], nNeg: Double = 0.1, kappa: Double = 1, nbIter: Int = 1000,
                    nbRestart: Int = 10, preMarkSum: Array[Double] = Array(), seed: Int = 0, verbose: Boolean = true): (List[Int], List[Int], Double) = {
        val expr = buildExprMap(m)
        val markSum = if (preMarkSum.length == m(0).length) preMarkSum else getMarkSum(m)

        val (initMarks, initObj) = getMarkers(m, initSamples.toList, expr, markSum, nNeg, kappa)
        var best = (initSamples, initMarks, initObj)

        if (verbose) {
            val i = ""
            val nb = "New best?"
            val c = "Nb. samples"
            val g = "Nb. markers"
            val o = "Obj. value"
            val t = "Exe. time [s]"
            println(f"| $i%-25s | $nb%-12s | $c%-12s | $g%-12s | $o%-12s | $t%-15s |")
            val i2 = "Initial solution"
            val nb2 = "    ****"
            val c2 = initSamples.size
            val g2 = initMarks.length
            val o2 = (initObj * 100).round / 100.0
            println(f"| $i2%-25s | $nb2%-12s | $c2%-12s | $g2%-12s | $o2%-12s | $i%-15s |")
        }

        for (nbr <- 1 to nbRestart) {
            if (verbose) {
                val i = "Local search restart " + nbr
                print(f"| $i%-25s |")
            }
            var bestRestart = (Set[Int](), List[Int](), 0.0)
            val t0 = System.currentTimeMillis

            var sp = Array.fill(m(0).length)(0.0)
            var sn = Array.fill(m(0).length)(0.0)
            var nn = Array.fill(m(0).length)(0)
            for (j <- m(0).indices) {
                for (i <- initSamples) {
                    if (m(i)(j) >= 0) {
                        sp(j) += m(i)(j)
                    } else {
                        sn(j) += m(i)(j)
                        nn(j) += 1
                    }
                }
            }
            var current = (initSamples, initMarks, initObj)

            val heurOrd = getSamObj(m, best._2).zipWithIndex.sortBy(-_._1).map(_._2)

            val rand = if (seed == 0) new Random() else new Random(seed * nbr)
            var temp = temperature(initObj, 0)

            for (k <- 0 until nbIter) {

                val (curOrdered, oocOrdered) = heurOrd.partition(current._1 contains _)

                val (samples, spNew, snNew, nnNew) = neighbour(m, curOrdered, oocOrdered, sp, sn, nn, rand)

                val candMark = nnNew.zipWithIndex.toList.filter { case (n, j) => n <= ceil(samples.size * nNeg) }
                        .map { case (n, j) => (j, n, -kappa * markSum(j) + (1 + kappa) * spNew(j) + snNew(j)) }
                        .filter(j => j._3 >= 0).sortBy(j => (j._2, -j._3))

                val (marks, obj) = selectCandMarkers(m, candMark, samples.size, nNeg)

                if (obj > bestRestart._3) {
                    bestRestart = (samples, marks, obj)

                    if (obj > best._3) {
                        best = (samples, marks, obj)
                    }
                }

                if (acceptance(current._3, obj, temp) >= rand.nextDouble()) {
                    current = (samples, marks, obj)
                    sp = spNew
                    sn = snNew
                    nn = nnNew
                }

                temp = temperature(initObj, k)
            }

            val (samples, marks, obj) = greedyImprove(m, bestRestart, markSum, nNeg, kappa)

            if (obj > bestRestart._3) {
                bestRestart = (samples, marks, obj)

                if (obj > best._3) {
                    best = (samples, marks, obj)
                }
            }

            val t1 = System.currentTimeMillis

            if (verbose) {
                val nb = if (bestRestart._3 >= best._3) "    ****" else ""
                val c = bestRestart._1.size
                val g = bestRestart._2.length
                val o = (bestRestart._3 * 100).round / 100.0
                val t = (t1 - t0).toDouble / 1000
                println(f" $nb%-12s | $c%-12s | $g%-12s | $o%-12s | $t%-15s |")
            }
        }

        return (best._1.toList.sorted, best._2, best._3)
    }

    /**
      * Sarcastically generates a new neighbour
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param curOrdered the samples currently selected, ordered decreasingly by sum of expression over initial markers
      * @param oocOrdered the samples currently not selected, ordered decreasingly by sum of expression over initial markers
      * @param sp sum of positive values in the current samples for each marker
      * @param sn sum of negative values in the current samples for each marker
      * @param nn number of negative values in the current samples for each marker
      * @param rand a pseudo-random number generator
      * @return a new assignment of cells, with updated sp, sn and nn
      */
    def neighbour(m: Array[Array[Double]], curOrdered: Array[Int], oocOrdered: Array[Int],
                  sp: Array[Double], sn: Array[Double], nn: Array[Int], rand: Random):
    (Set[Int], Array[Double], Array[Double], Array[Int]) = {

        if (rand.nextInt(m.length) < oocOrdered.length) { // Add sample to current
            val idx = floor(oocOrdered.length * pow(rand.nextDouble(), 2)).toInt
            val sam = oocOrdered(idx)

            val spNew = (sp.indices).toArray.map(j => {
                val a = if (m(sam)(j) >= 0) m(sam)(j) else 0.0
                sp(j) + a
            })
            val snNew = (sn.indices).toArray.map(j => {
                val a = if (m(sam)(j) < 0) m(sam)(j) else 0.0
                sn(j) + a
            })
            val nnNew = (nn.indices).toArray.map(j => {
                val a = if (m(sam)(j) < 0) 1 else 0
                nn(j) + a
            })

            return (curOrdered.toSet + sam, spNew, snNew, nnNew)

        } else { // Remove sample from current
            val idx = (curOrdered.length - 1) - floor(curOrdered.length * pow(rand.nextDouble(), 2)).toInt
            val sam = curOrdered(idx)

            val spNew = (sp.indices).toArray.map(j => {
                val a = if (m(sam)(j) >= 0) m(sam)(j) else 0.0
                sp(j) - a
            })
            val snNew = (sn.indices).toArray.map(j => {
                val a = if (m(sam)(j) < 0) m(sam)(j) else 0.0
                sn(j) - a
            })
            val nnNew = (nn.indices).toArray.map(j => {
                val a = if (m(sam)(j) < 0) 1 else 0
                nn(j) - a
            })

            return (curOrdered.toSet - sam, spNew, snNew, nnNew)
        }
    }

    /**
      * Temperature for the simulated annealing
      * @param t0 the initial temperature
      * @param k the iteration number
      * @return the current temperature
      */
    def temperature(t0: Double, k: Int): Double = {
        t0 * pow(0.99, k)
    }

    /**
      * Computes the probability to accept the neighbour
      * @param obj the current objective value
      * @param objNew the objective value of the neighbour
      * @param temp the current temperature
      * @return a probability
      */
    def acceptance(obj: Double, objNew: Double, temp: Double): Double = {
        if (objNew >= obj) {
            return 1.0
        } else {
            return (exp(-(obj - objNew) / temp)) max 0.01
        }
    }

    /**
      * Greedily searches if a better solution can be found by adding out-of-cluster samples, in order of their sum
      * of expression over the selected markers
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param init the current solution
      * @param markSum the sum of positive values for each marker
      * @param nNeg the maximum percentage of negative values allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return a (possibly improved) solution
      */
    def greedyImprove(m: Array[Array[Double]], init: (Set[Int], List[Int], Double), markSum: Array[Double],
                      nNeg: Double, kappa: Double): (Set[Int], List[Int], Double) = {
        var current = init._1
        var best = init

        var sp = Array.fill(m(0).length)(0.0)
        var sn = Array.fill(m(0).length)(0.0)
        var nn = Array.fill(m(0).length)(0)
        for (j <- m(0).indices) {
            for (i <- init._1) {
                if (m(i)(j) >= 0) {
                    sp(j) += m(i)(j)
                } else {
                    sn(j) += m(i)(j)
                    nn(j) += 1
                }
            }
        }

        val samObj = getSamObj(m, init._2).zipWithIndex
        val (inclSamUnsorted, oocSamUnsorted) = samObj.partition(init._1 contains _._2)
        val inclSam = inclSamUnsorted.sortBy(_._1)
        val oocSam = oocSamUnsorted.sortBy(-_._1)

        for (i <- 0 until (oocSam.length min inclSam.length)) {
            if (oocSam(i)._1 > inclSam(i)._1) {
                current = current + oocSam(i)._2 - inclSam(i)._2

                sp = (sp.indices).toArray.map(j => {
                    val a = if (m(oocSam(i)._2)(j) >= 0) m(oocSam(i)._2)(j) else 0.0
                    val b = if (m(inclSam(i)._2)(j) >= 0) m(inclSam(i)._2)(j) else 0.0
                    sp(j) + a - b
                })
                sn = (sn.indices).toArray.map(j => {
                    val a = if (m(oocSam(i)._2)(j) < 0) m(oocSam(i)._2)(j) else 0.0
                    val b = if (m(inclSam(i)._2)(j) < 0) m(inclSam(i)._2)(j) else 0.0
                    sn(j) + a - b
                })
                nn = (nn.indices).toArray.map(j => {
                    val a = if (m(oocSam(i)._2)(j) < 0) 1 else 0
                    val b = if (m(inclSam(i)._2)(j) < 0) 1 else 0
                    nn(j) + a - b
                })
            } else {
                current = current + oocSam(i)._2

                sp = (sp.indices).toArray.map(j => {
                    val a = if (m(oocSam(i)._2)(j) >= 0) m(oocSam(i)._2)(j) else 0.0
                    sp(j) + a
                })
                sn = (sn.indices).toArray.map(j => {
                    val a = if (m(oocSam(i)._2)(j) < 0) m(oocSam(i)._2)(j) else 0.0
                    sn(j) + a
                })
                nn = (nn.indices).toArray.map(j => {
                    val a = if (m(oocSam(i)._2)(j) < 0) 1 else 0
                    nn(j) + a
                })
            }

            val candMark = nn.zipWithIndex.toList.filter { case (n, j) => n <= ceil(current.size * nNeg) }
                    .map { case (n, j) => (j, n, -kappa * markSum(j) + (1 + kappa) * sp(j) + sn(j)) }
                    .filter(j => j._3 >= 0).sortBy(j => (j._2, -j._3))

            val (marks, obj) = selectCandMarkers(m, candMark, current.size, nNeg)

            if (obj > best._3) {
                best = (current, marks, obj)
            }
        }

        return best
    }
}
