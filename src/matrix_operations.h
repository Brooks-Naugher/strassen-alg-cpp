#include <vector>

/**
 * @brief Adds two matrix vectors (1d).
 * 
 * Adds two, one dimentional, same-size vectors, together,
 * handling it like matrix addition.
 * 
 * @param matrixA
 * @param matrixB
 * @return resultMatrix - sum of a & b
 * 
 */
std::vector<int> addMatrix(const std::vector<int>& matrixA, const std::vector<int>& matrixB){
    std::vector<int> resultMatrix(matrixA.size(), 0);

    for(std::size_t i = 0; i < matrixA.size(); i++){
        resultMatrix[i] = matrixA[i] + matrixB[i];
    }

    return resultMatrix;
}

/**
 * @brief Subtracts two matrix vectors (1d).
 * 
 * Subtracts two, one dimentional, same-size vectors, 
 * handling it like matrix subtraction.
 * 
 * @param matrixA
 * @param matrixB
 * @return resultMatrix - subtraction of a & b
 * 
 */
std::vector<int> subtractMatrix(const std::vector<int>& matrixA, const std::vector<int>& matrixB){
    std::vector<int> resultMatrix(matrixA.size(), 0);

    for(std::size_t i = 0; i < matrixA.size(); i++){
        resultMatrix[i] = matrixA[i] - matrixB[i];
    }

    return resultMatrix;
}

