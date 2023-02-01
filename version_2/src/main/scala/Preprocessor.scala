/*
 * Author: Alexander Gerniers
 */
object Preprocessor {

    /**
      * Filters out the markers that are expressed in only 0 or 1 samples,
      * and those expressed in more than a certain proportion thresh of the samples
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param markNames the list of names corresponding to each marker
      * @param thresh in [0, 1]
      * @return an expression matrix and marker name array
      */
    def expressionFilter(m: Array[Array[Double]], markNames: Array[String], thresh: Double): (Array[Array[Double]], Array[String]) = {
        val nSam = m.length
        val mt = m.transpose
        val nb = mt.map(l => m.length - l.groupBy(x => x.toInt).getOrElse(-1, Array()).length).zipWithIndex
        val idx =  nb.filter(_._1 <= thresh * nSam).filterNot(_._1 <= 1).map(_._2)
        return ((idx map mt).transpose, idx map markNames)
    }

}
