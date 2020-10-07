import scala.collection.mutable.ListBuffer

import Objective.{getMarkers, buildExprMap, getMarkSum, getCoExprMarks}

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
      * @param deltaHeuristic the number of extra top-solutions to consider at each level
      * @param maxNbSam the maximum number of samples in the bicluster
      * @param minCoExpMark a minimum number of markers that must be expressed in every row of the bicluster
      * @param excl a list of markers that need to be excluded from the search
      * @return - an assignment of samples
      *         - an assignment of markers
      *         - the corresponding objective value
      */
    def findCluster(m: Array[Array[Double]], nNeg: Double = 0.1, kappa: Double = 1, nHeuristic: Int = 20, deltaHeuristic: Int = 10,
                    maxNbSam: Int = Int.MaxValue, minCoExpMark: Int = 0, excl: List[Int] = List()): (List[Int], List[Int], Double) = {
        val expr = buildExprMap(m)
        val markSum = getMarkSum(m)
        val nSam = m.length

        val lb = ListBuffer[(List[Int], List[Int], Double)]()
        val pairs = (0 until nSam).combinations(2).toList

        println("Current search level: 2")
        println("\t " + pairs.length + " pairs to evaluate")

        // Evaluation of all the pairs of cells
        for (p <- pairs) {
            if (excl.length == 0 || !excl.exists(p contains _)) {
                val (markers, obj) = getMarkers(m, p.toList, expr, markSum, nNeg, kappa)
                lb += ((p.toList, markers, obj))
            }
        }

        val l = lb.sortBy(-_._3).toList.take(nHeuristic)
        var best = l(0)
        var prevBest = l.map(x =>  (x._1, getCoExprMarks(x._1, expr)))

        println("\t NEW BEST: obj. val.: " + best._3 + "\t samples (idx. from 1): " + best._1.map(_ + 1).mkString(", ") + "\t nb. markers: " + best._2.length)

        var lvl = 3
        var finished = false
        while (!finished && lvl <= maxNbSam) {
            // Extend solutions in prevBest by adding one cell
            var candidates = prevBest.flatMap(x => (0 until nSam).map(y => ((y :: x._1).distinct, (x._2.toSet intersect expr(y).toSet).toList)))

            // Filter out solutions containing excluded cells
            if (excl.length > 0) candidates = candidates.filterNot(x => excl.exists(x._1 contains _))

            // Filter out solutions with smaller number of samples and without enough markers expressed in all samples
            candidates = candidates.filterNot(_._1.size < lvl)
                                   .map(x => (x._1.sorted, x._2)).distinct
                                   .filter(x => x._2.length >= minCoExpMark)

            println("Current search level: " + lvl)
            println("\t " + candidates.length + " solutions to evaluate")

            if (!candidates.isEmpty) {
                var res = ListBuffer[(List[Int], List[Int], Double, List[Int])]()

                // Evaluation of the solutions
                for ((samples, coExp) <- candidates) {
                    val (markers, obj) = getMarkers(m, samples, expr, markSum, nNeg, kappa)
                    res += ((samples, markers, obj, coExp))
                }

                if (res.length > 0) {
                    res = res.sortBy(-_._3)
                    prevBest = res.toList.take(nHeuristic + (lvl - 2) * deltaHeuristic).map(x => (x._1, x._4))
                    if (res(0)._3 > best._3) {
                        best = (res(0)._1, res(0)._2, res(0)._3)

                        println("\t NEW BEST: obj. val.: " + best._3 + "\t samples (idx. from 1): " + best._1.map(_ + 1).mkString(", ") + "\t nb. markers: " + best._2.length)
                    }
                    lvl += 1
                } else {
                    finished = true
                }
            } else {
                finished = true
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
