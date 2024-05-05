import java.io.File
import kotlin.math.abs

data class CalculationResult(
    val result: DoubleArray,
    val iterations: Int
) {
    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false

        other as CalculationResult

        if (!result.contentEquals(other.result)) return false
        if (iterations != other.iterations) return false

        return true
    }

    override fun hashCode(): Int {
        var result1 = result.contentHashCode()
        result1 = 31 * result1 + iterations
        return result1
    }
}

const val MAX_ITERATIONS = 1000000

fun gaussSeidel(
    coefficients: Array<DoubleArray>,
    rhs: DoubleArray,
    tolerance: Double
): CalculationResult? {
    val numberOfEquations = coefficients.size
    val result = DoubleArray(numberOfEquations)
    val maxIterations = MAX_ITERATIONS
    val previousX = DoubleArray(numberOfEquations)
    for (i in 0..<maxIterations) {
        result.forEachIndexed { index, _ ->
            previousX[index] = result[index]
        }
        for (j in 0..<numberOfEquations) {
            var sum = 0.0
            for (k in 0..<numberOfEquations) {
                if (k != j) {
                    sum += coefficients[j][k] * result[k]
                }
            }
            result[j] = (rhs[j] - sum) / coefficients[j][j]
        }
        var diff1norm = 0.0
        var oldNorm = 0.0
        for (j in 0..<numberOfEquations) {
            diff1norm += abs(result[j] - previousX[j])
            oldNorm += abs(previousX[j])
        }
        if (oldNorm == 0.0) {
            oldNorm = 1.0
        }
        val norm = diff1norm / oldNorm
        if (norm < tolerance && i != 0) {
            return CalculationResult(result, i + 1)
        }
    }
    return null
}


fun getEpsilons(
    coefficients: Array<DoubleArray>,
    solution: DoubleArray,
    rhs: DoubleArray,
): DoubleArray {
    val numberOfEquations = coefficients.size
    val result = DoubleArray(numberOfEquations)
    for (rowIndex in coefficients.indices) {
        var sum = 0.0
        for (columnIndex in coefficients[rowIndex].indices) {
            sum += coefficients[rowIndex][columnIndex] * solution[columnIndex]
        }
        result[rowIndex] = rhs[rowIndex] - sum
    }
    return result
}

fun matrixFromFile(filename: String): Array<DoubleArray> {
    val file = File(filename)
    val listMatrix = file.useLines { it.toList() }.map { it.split(" ").map { it.toDouble() } }
    val result = listMatrix.mapIndexed { idx, element ->
        val lol = DoubleArray(element.size)
        for (listElementIdx in element.indices) {
            lol[listElementIdx] = element[listElementIdx]
        }
        lol
    }.toTypedArray()
    return result
}

fun extractColumn(columnIndex: Int, matrix: Array<DoubleArray>): DoubleArray {
    val size = matrix.size
    val result = DoubleArray(size)
    for (i in 0..<size) {
        result[i] = matrix[i][columnIndex]
    }
    return result
}


fun Array<DoubleArray>.removeLastColumn() {
    for (rowIdx in this.indices) {
        this[rowIdx] = this[rowIdx].slice(0..<this[rowIdx].lastIndex).toDoubleArray()
    }
}


fun main() {
    val matrix = matrixFromFile("matrix.txt")
    val rhs = extractColumn(matrix.size, matrix)
    matrix.removeLastColumn()

    val result2 = gaussSeidel(matrix, rhs, 1E-14)
    if (result2 == null) {
        println("Результат рассчитать не удалось!")
        return
    }
    println("Результат полученный за ${result2.iterations} итераций: ${result2.result.toList()}")
    println("Соответствующая невязка: ${getEpsilons(matrix, result2.result, rhs).toList()}")
}