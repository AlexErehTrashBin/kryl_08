import java.io.File
import java.util.*
import kotlin.collections.ArrayList
import kotlin.math.abs
import kotlin.math.max
import kotlin.system.exitProcess

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

fun isDiagonalPriority(a: Array<DoubleArray>, n: Int): Boolean {
    var flag = true

    for (i in 0..<n) {
        var s = 0.0
        for (j in 0..<n) {
            s += abs(a[i][j])
        }
        s -= abs(a[i][i])
        if (s > abs(a[i][i])) {
            flag = false
        }
    }

    return flag
}

//    пытаемся сделать диагональное преобладание
fun doDiagonalPriority(n: Int, a: Array<DoubleArray>, b: DoubleArray): Array<DoubleArray> {
    val options: MutableList<SortedSet<Int>> = ArrayList()
    for (i in 0..<n) {
        options.add(TreeSet())
    }

    var k = 0
    for (i in 0..<n) {
        var s = 0.0
        for (j in 0..<n) {
            s += abs(a[i][j])
        }
        for (j in 0..<n) {
            if (abs(a[i][j]) >= (s - abs(a[i][j]))) {
                val set = options[j]
                set.add(i)
            }
            if (abs(a[i][j]) > (s - abs(a[i][j]))) {
                k++
            }
        }
    }

    if (k == 0) {
        throw RuntimeException("Невозможно преобразовать к диагональному виду ;(")
    }
    checkPotentialDiagonalPriority(options, n)

    var maxLen = 20
    while (maxLen > 2) {
        var curMaxLen = 0
        for (i in 0..<n) {
            for (j in 0..<n) {
                if (i != j && options[j].size == 1) {
                    options[i].removeAll(options[j])
                    curMaxLen = max(curMaxLen.toDouble(), options[i].size.toDouble()).toInt()
                }
            }
        }
        maxLen = curMaxLen
    }

    checkPotentialDiagonalPriority(options, n)

    //        получаем необходимую последовательность строк
    for (i in 0..<n) {
        if (options[i].size == 2) {
            options[i].remove(options[i].last())
        }
        for (j in i + 1..<n) {
            options[j].removeAll(options[i])
        }
    }

    checkPotentialDiagonalPriority(options, n)

    val newArray = Array(n) { DoubleArray(n + 1) }

    for (i in 0..<n) {
        for (j in 0..<n) {
            newArray[i][j] = a[options[i].first()][j]
        }
        newArray[i][n] = b[options[i].first()]
    }

    return newArray
}

private fun checkPotentialDiagonalPriority(options: List<SortedSet<Int>>, n: Int) {
    for (i in 0..<n) {
        if (options[i].isEmpty()) {
            throw RuntimeException("Невозможно преобразовать к диагональному виду :(")
        }
    }
}

fun gaussSeidel(
    coefficients: Array<DoubleArray>,
    rhs: DoubleArray,
    tolerance: Double
): CalculationResult? {
    try {
        if (!isDiagonalPriority(coefficients, coefficients.size)) {
            val newA = doDiagonalPriority(coefficients.size, coefficients, rhs)
            for (i in coefficients.indices) {
                System.arraycopy(newA[i], 0, coefficients[i], 0, coefficients.size)
                rhs[i] = newA[i][coefficients.size]
            }
        }
    } catch (e: RuntimeException) {
        println(e.message)
        exitProcess(1)
    }

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