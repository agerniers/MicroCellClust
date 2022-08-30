import scala.collection.mutable.{ListBuffer, Map => MutMap}
import scala.math.ceil

/*
 * Author: Alexander Gerniers (UCLouvain)
 */
object Objective {

    /**
      * Get the assignment of markers given the assignment of samples
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param samples an assignment of samples
      * @param expr map containing for each row the expressed columns
      * @param markSum the cumulative positive expression of the markers
      * @param nNeg the maximum percentage of negative values allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return - the assignment of markers respecting the maximal percentage of negative values using a greedy approach
      *           (cols ranked per number of negative values, then per obj value)
      *         - the corresponding objective value
      */
    def getMarkers(m: Array[Array[Double]], samples: List[Int], expr: Map[Int, List[Int]], markSum: Array[Double],
                   nNeg: Double = 0.1, kappa: Double = 1): (List[Int], Double) = {
        val markExpr = samples.flatMap(expr(_))

        val candidateMarks = markExpr.groupBy(identity)
                .map(j => (j._1, samples.length - j._2.size, getMarkObj(m, samples, j._1, markSum, kappa)))
                .toList.filter(j => j._3 >= 0 && j._2 <= ceil(samples.length * nNeg))
                .sortBy(j => (j._2, -j._3))

        return selectCandMarkers(m, candidateMarks, samples.size, nNeg)
    }

    /**
      * Get the assignment of markers given the assignment of samples
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param nbSam the number of samples in the assignment
      * @param newSample the sample that was added
      * @param markSum the cumulative positive expression of the markers
      * @param prevExpSums a marker -> (sum pos, sum neg, nb neg) map of the previous sample assignment
      * @param nNeg the maximum percentage of negative values allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return - the assignment of markers respecting the maximal percentage of negative values using a greedy approach
      *           (cols ranked per number of negative values, then per obj value)
      *         - the corresponding objective value
      */
    def getMarkersFromPrev(m: Array[Array[Double]], nbSam: Int, newSample: Int, markSum: Array[Double],
                           prevExpSums: List[(Int, Double, Double, Int)], nNeg: Double = 0.1, kappa: Double = 1): (List[Int], Double) = {
        val maxNeg = ceil(nbSam * nNeg)

        val candidateMarks = prevExpSums.map{case (j, sp, sn, nn) => getMarkObjFromPrev(m, newSample, j, markSum, sp, sn, nn, kappa)}
                .filter(j => j._3 >= 0 && j._2 <= maxNeg)
                .sortBy(j => (j._2, -j._3))

        return selectCandMarkers(m, candidateMarks, nbSam, nNeg)
    }

    /**
      * Selects the assignment of markers from a list of candidates for a certain assignment of samples
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param candidateMarks a list of markers of format (idx_marker, nb_negatives, obj_val)
      * @param nbSam the number of samples in the assignment
      * @param nNeg the maximum percentage of negatives allowed inside the cluster
      * @return a list of markers and corresponding objective value
      */
    def selectCandMarkers(m: Array[Array[Double]], candidateMarks: List[(Int, Int, Double)], nbSam: Int, nNeg: Double): (List[Int], Double) = {
        val maxNeg = ceil(nbSam * nNeg)

        val markers = new ListBuffer[Int]
        var obj = 0.0
        var finished = false
        var n = 0.0
        var j = 0
        while (!finished && j < candidateMarks.length) {
            val mark = candidateMarks(j)
            if ((n + mark._2) / (nbSam * (markers.length + 1)) <= nNeg) {
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
      * @param markSum the cumulative positive expression of the markers
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return the objective value
      */
    def getMarkObj(m: Array[Array[Double]], samples: List[Int], marker: Int, markSum: Array[Double], kappa: Double): Double = {
        var obj = - kappa * markSum(marker)

        for (i <- samples) {
            if (m(i)(marker) >= 0) obj += (1 + kappa) * m(i)(marker)
            else obj += m(i)(marker)
        }

        return obj
    }

    /**
      * Gives the contribution of a marker to the objective when adding one sample to the previous assignment
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param newSample the sample that was added
      * @param marker a marker
      * @param kappa a weighting constant for the out-of-cluster expression
      * @param markSum the cumulative positive expression of the markers
      * @param prevPosSum sum of positive expression of the marker in the previous sample assignment
      * @param prevNegSum sum of negative expression of the marker in the previous sample assignment
      * @param prevNegNb number of negative expression of the marker in the previous sample assignment
      * @return the objective value and number of negatives : (marker, nn, obj)
      */
    def getMarkObjFromPrev(m: Array[Array[Double]], newSample: Int, marker: Int, markSum: Array[Double],
                           prevPosSum: Double, prevNegSum: Double, prevNegNb: Int, kappa: Double): (Int, Int, Double) = {
        var obj = - kappa * markSum(marker) + (1 + kappa) * prevPosSum + prevNegSum
        var nn = prevNegNb

        if (m(newSample)(marker) >= 0) {
            obj += (1 + kappa) * m(newSample)(marker)
        } else {
            obj += m(newSample)(marker)
            nn += 1
        }

        return (marker, nn, obj)
    }

    /**
      * Gives the sum of positive, sum and number of negative expression of markers for the specified sample assignment
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param samples an assignment of samples
      * @return A (marker, sum pos, sum neg, nb neg) list
      */
    def getExpSums(m: Array[Array[Double]], samples: List[Int], nNeg: Double): List[(Int, Double, Double, Int)] = {
        val expSums = ListBuffer[(Int, Double, Double, Int)]()
        val maxNegNext = ceil((samples.length + 1) * nNeg)

        for (j <- m(0).indices) {
            var sp = 0.0
            var sn = 0.0
            var nn = 0
            for (i <- samples) {
                if (m(i)(j) >= 0) {
                    sp += m(i)(j)
                } else {
                    sn += m(i)(j)
                    nn += 1
                }
            }
            if (nn <= maxNegNext) {
                expSums += ((j, sp, sn, nn))
            }
        }

        return expSums.toList
    }

    /**
      * Builds a mapping containing, for each sample of m, the markers that are expressed
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @return a sample -> list of markers map
      */
    def buildExprMap(m: Array[Array[Double]]): Map[Int, List[Int]] = {
        val expr = MutMap[Int, List[Int]]()

        for (i <- m.indices) {
            expr += (i -> m(i).zipWithIndex.filter(_._1 >= 0).map(_._2).toList)
        }

        return expr.toMap
    }

    /**
      * Get the cumulative positive expression of each marker
      * @param m an expression matrix
      * @return an array with the sum of positive expressions of each marker
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
      * Computes the percentage of the given markers expressed by the sample
      * @param sample a sample
      * @param markers a list of markers
      * @param expr a sample -> markers map
      * @return the percentage
      */
    def getSupport(sample: Int, markers: List[Int], expr: Map[Int, List[Int]]): Double = {
        val inter = markers intersect expr(sample)
        return inter.length.toDouble / markers.length
    }

    /**
      * Computes the sum of expression over given markers for each sample
      * @param m an expression matrix
      * @param markers a list of markers
      * @return an array with the sum of maker expression for each sample
      */
    def getSamObj(m: Array[Array[Double]], markers: List[Int]): Array[Double] = {
        return m.indices.toArray.map(i => markers.map(j => m(i)(j)).sum)
    }
}
