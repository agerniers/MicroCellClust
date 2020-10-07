import scala.collection.mutable.ListBuffer
import scala.io.Source

/*
 * Author: Alexander Gerniers
 */
object DataReader {

    /**
      * Reads the expression values contained in a text file
      * @param file the path to the file to be converted
      * @param transpose whether to transpose the matrix (recall: samples (cells) on rows and markers (genes) on columns)
      * @param names whether the first line and column in the file contains names
      * @return a 2D array containing the expression values
      */
    def readMatrix(file: String, transpose: Boolean = false, names: Boolean = false): Array[Array[Double]] = {

        var m = new ListBuffer[Array[Double]]

        if (names) {
            for (line <- Source.fromFile(file).getLines.drop(1)) {
                m += line.split("\\s+|,").filterNot(_ == "").drop(1).map(_.toDouble)
            }
        } else {
            for (line <- Source.fromFile(file).getLines) {
                m += line.split("\\s+|,").filterNot(_ == "").map(_.toDouble)
            }
        }

        return if (transpose) m.toArray.transpose else m.toArray
    }

    /**
      * Reads a file containing the names of the markers (genes) on each line
      * @param file the filepath
      * @return an array of strings
      */
    def readMarkNames(file: String): Array[String] = Source.fromFile(file).getLines.toArray

}
