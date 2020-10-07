import scala.collection.mutable.{ListBuffer, Map => MutMap}
import scala.math.ceil

/*
 * Author: Alexander Gerniers
 */
object Objective {

    /**
      * Gives the contribution of a marker to the objective under the given assignment of samples
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param samples an assignment of samples
      * @param expr map containing for each row the expressed columns
      * @param markSum the cumulative positive expression of the markers
      * @param nNeg the maximum percentage of -1 allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return - the assignment of markers respecting the maximal percentage of -1 using a greedy approach
      *           (cols ranked per number of -1, then per obj value)
      *         - the corresponding objective value
      */
    def getMarkers(m: Array[Array[Double]], samples: List[Int], expr: Map[Int, List[Int]], markSum: Array[Double],
                   nNeg: Double = 0.1, kappa: Double = 1): (List[Int], Double) = {
        val markExpr = new ListBuffer[Int]()
        samples.foreach(markExpr ++= expr(_))

        val candidateMarks = markExpr.toList.groupBy(identity)
                .map(j => (j._1, samples.length - j._2.size, getMarkObj(m, samples, j._1, kappa, markSum)))
                .toList.filter(j => j._3 >= 0 && j._2 <= ceil(samples.length * nNeg)).sortBy(j => (j._2, -j._3))

        val markers = new ListBuffer[Int]
        var obj = 0.0
        var finished = false
        var n = 0.0
        var j = 0
        while (!finished && j < candidateMarks.length) {
            val mark = candidateMarks(j)
            if ((n + mark._2) / (samples.length * (markers.length + 1)) <= nNeg) {
                markers += mark._1
                n += mark._2
                obj += mark._3
                j += 1
            } else finished = true
        }

        return (markers.toList, obj)
    }

    /**
      * Gives the contribution of a marker to the objective under the given assignment of samples
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param samples an assignment of samples
      * @param marker a marker
      * @param kappa a weighting constant for the out-of-cluster expression
      * @param markSum the cumulative positive expression of the markers
      * @return the objective value
      */
    def getMarkObj(m: Array[Array[Double]], samples: List[Int], marker: Int, kappa: Double = 1, markSum: Array[Double]): Double = {
        var obj = - kappa * markSum(marker)

        for (i <- samples) {
            if (m(i)(marker) >= 0) obj += (1 + kappa) * m(i)(marker)
            else obj += m(i)(marker)
        }

        return obj
    }

    /**
      * Builds a mapping containing, for each sample of m, the markers that are expressed
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @return a sample -> marker map
      */
    def buildExprMap(m: Array[Array[Double]]): Map[Int, List[Int]] = {
        val expr = MutMap[Int, List[Int]]()

        for (i <- m.indices) {
            expr += (i -> m(i).zipWithIndex.filter(_._1 >= 0).map(_._2).toList)
        }

        return expr.toMap
    }

    /**
      * Get the cumulative positive expression of each column
      * @param m an expression matrix
      * @return an array with the sum of positive expressions of each column
      */
    def getMarkSum(m: Array[Array[Double]]): Array[Double] = {
        val markSum = Array.fill(m(0).length){0.0}

        for (j <- m(0).indices) {
            for (i <- m.indices) {
                if (m(i)(j) >= 0) markSum(j) += m(i)(j)
            }
        }

        return markSum
    }

    /**
      * Get the markers that are expressed for all the samples in the given assignment
      * @param samples an assignment of samples
      * @param expr a sample -> marker map
      * @return an assignment of markers
      */
    def getCoExprMarks(samples: List[Int], expr: Map[Int, List[Int]]): List[Int] = {
        var interMarks = expr(samples(0)).toSet

        for (i <- 1 until samples.length) {
            interMarks = interMarks intersect expr(samples(i)).toSet
        }

        return interMarks.toList.sorted
    }
}
