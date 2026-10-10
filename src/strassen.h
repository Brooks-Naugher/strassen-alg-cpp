#include <vector>
#include <cmath>
#include <type_traits>
#include "matrix_operations.h"

//https://en.wikipedia.org/wiki/Strassen_algorithm


/**
 * @brief Checks if a matrix can be multiplied with Strassen Algorithm.
 * 
 * Checks if log2 of square root of number of elements in
 * the matrix is a whole number, if it isn't, then the width/length 
 * is not a power of 2 and/or square.
 * 
 * @param matrix1D 
 * @return bool
 * 
 */
bool isStrassenMultipliable1d(std::vector<int>& matrix1D){
    
    size_t size = matrix1D.size();
    
    if(size == 0){
        return false;
    }

    if(std::fmod(std::log2(sqrt(size)), 1.0) != 0 ){
        return false;
    }

    return true;
}

/**
 * @brief Checks if two matrices can be multiplied with Strassen Algorithm with eachother.
 * 
 * Checks if size of matrices are same and that both matrices meet requirements of being
 * Strassen multipliable.
 * 
 * @param matrixA 
 * @param matrixB 
 * @return bool
 * 
 */
bool areBothStrassenMultipliable(std::vector<int>& matrixA, std::vector<int>& matrixB){
    size_t matrixASize = matrixA.size();
    size_t matrixBSize = matrixB.size();

    bool areSameSize = matrixASize == matrixBSize;
    if(!areSameSize){
        return false;
    }

    bool areStrassenMultipliable = isStrassenMultipliable1d(matrixA) && isStrassenMultipliable1d(matrixB);
    
    if(!areStrassenMultipliable){
        return false; 
    }

    return true;

}
/**
 * @brief Multiplies two corner ordered matrices using Strassen Algorithm.
 * 
 * Multiplies two same-sized, Strassen multipliable
 * matrices. The matrices are in "Corner Order".
 * Corner Order Matrix  looks like this: If
        m = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
        Then each corner: 
        m11 = 1,2,3,4 | 
        m12 = 4,6,7,8 | 
        m21 = 9,10,11,12 | 
        m22 = 13,14,15,16
    @param matrixA
    @param matrixB
    @return matrixC - multiplied result of matrixA and matrixB

 */
std::vector<int> strassenMultiplyCornerOrder(const std::vector<int>& matrixA, const std::vector<int>& matrixB){

    /*
            Corner Order Matrix  looks like this
            m = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}
            m11 = 1,2,3,4
            m12 = 4,6,7,8
            m21 = 9,10,11,12
            m22 = 13,14,15,16
    */

    std::size_t numElements = matrixA.size();
    std::vector<int> strassenResult;

    if(numElements == 1){
        strassenResult = {matrixA.at(0)*matrixB.at(0)};
    }
    else if(numElements != 4){

        std::size_t numCornerElements = numElements / 4;
        
        std::vector<int> a11, a12, a21, a22;
        std::vector<int> b11, b12, b21, b22;
        std::vector<int> m1, m2, m3, m4, m5, m6, m7;
        std::vector<int> c11, c12, c21, c22;

        for(std::size_t i = 0; i < numElements; i++){
            if(i < numCornerElements){
                a11.push_back(matrixA.at(i));
                b11.push_back(matrixB.at(i));
            } 
            else if(i < numCornerElements * 2){
                a12.push_back(matrixA.at(i));
                b12.push_back(matrixB.at(i));
            }
            else if(i < numCornerElements * 3){
                a21.push_back(matrixA.at(i));
                b21.push_back(matrixB.at(i));
            }
            else{
                a22.push_back(matrixA.at(i));
                b22.push_back(matrixB.at(i)); 
            }
        }

    /*
        m1 = (a11 + a22) * (b11 + b22);
        m2 = (a21 + a22) * b11;
        m3 = a11 * (b12-b22);
        m4 = a22 * (b21 - b11);
        m5 = (a11 + a12) * b22;
        m6 = (a21 - a11) * (b11 + b12);
        m7 = (a12 - a22) * (b21 + b22);
    */
        m1 = strassenMultiplyCornerOrder(addMatrix(a11, a22), addMatrix(b11, b22));
        m2 = strassenMultiplyCornerOrder(addMatrix(a21, a22), b11);
        m3 = strassenMultiplyCornerOrder(a11, subtractMatrix(b12, b22));
        m4 = strassenMultiplyCornerOrder(a22, subtractMatrix(b21, b11));
        m5 = strassenMultiplyCornerOrder(addMatrix(a11, a12), b22);
        m6 = strassenMultiplyCornerOrder(subtractMatrix(a21,a11), addMatrix(b11, b12));
        m7 = strassenMultiplyCornerOrder(subtractMatrix(a12, a22), addMatrix(b21, b22));

    /*
        c11 = m1 + m4 - m5 + m7;
        c12 = m3 + m5;
        c21 = m2 + m4;
        c22 = m1 - m2 + m3 + m6;
    */
        c11 = subtractMatrix(addMatrix(addMatrix(m1, m4), m7), m5);
        c12 = addMatrix(m3, m5);
        c21 = addMatrix(m2, m4);
        c22 = subtractMatrix(addMatrix(addMatrix(m1,m3),m6), m2);

        strassenResult.insert(strassenResult.end(), c11.begin(), c11.end());
        strassenResult.insert(strassenResult.end(), c12.begin(), c12.end());
        strassenResult.insert(strassenResult.end(), c21.begin(), c21.end());
        strassenResult.insert(strassenResult.end(), c22.begin(), c22.end());
    }
    else{

        int a11, a12, a21, a22;
        int b11, b12, b21, b22;
        int m1, m2, m3, m4, m5, m6, m7;
        int c11, c12, c21, c22;

        a11 = matrixA.at(0);
        b11 = matrixB.at(0);
        a12 = matrixA.at(1);
        b12 = matrixB.at(1);
        a21 = matrixA.at(2);
        b21 = matrixB.at(2);
        a22 = matrixA.at(3);
        b22 = matrixB.at(3);

        m1 = (a11 + a22) * (b11 + b22);
        m2 = (a21 + a22) * b11;
        m3 = a11 * (b12-b22);
        m4 = a22 * (b21 - b11);
        m5 = (a11 + a12) * b22;
        m6 = (a21 - a11) * (b11 + b12);
        m7 = (a12 - a22) * (b21 + b22);    
        
        c11 = m1 + m4 - m5 + m7;
        c12 = m3 + m5;
        c21 = m2 + m4;
        c22 = m1 - m2 + m3 + m6;

        strassenResult = {c11, c12, c21, c22};
    }

    return strassenResult;
        
}


/**
 * @brief Multiplies two row ordered matrices using Strassen Algorithm.
 * 
 * Multiplies two same-sized, Strassen multipliable
 * matrices. The matrices are in "Row Order".
 * Row Order Matrix  looks like this: If
        m = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
        Each Row: 
        m1 = 1,2,3,4 | 
        m2 = 4,6,7,8 | 
        m3 = 9,10,11,12 | 
        m4 = 13,14,15,16
    @param matrixA
    @param matrixB
    @return matrixC - multiplied result of matrixA and matrixB

 */

std::vector<int> strassenMultiplyRowOrder(const std::vector<int>& matrixA, const std::vector<int>& matrixB){
    

    std::size_t numElements = matrixA.size();
    std::vector<int> strassenResult(numElements, 0);

    if(numElements == 1){
        strassenResult = {matrixA.at(0)*matrixB.at(0)};
    }
    else if(numElements != 4){

        std::size_t numCornerElements = numElements / 4;
        std::vector<int> a11, a12, a21, a22;
        std::vector<int> b11, b12, b21, b22;
        std::vector<int> m1, m2, m3, m4, m5, m6, m7;
        std::vector<int> c11, c12, c21, c22;

        std::size_t rowNum = 1;
        std::size_t counter = 0;
        
        for(std::size_t i = 0; i < numElements; i++){
            if(i < numElements/2){
                bool isOfElement11 = i < (numCornerElements * rowNum) - (numCornerElements/2);
                if(isOfElement11){
                    a11.push_back(matrixA.at(i));
                    b11.push_back(matrixB.at(i));
                } 
                else {
                    a12.push_back(matrixA.at(i));
                    b12.push_back(matrixB.at(i));                }
            }
            else{
                bool isOfElement21 = i < (numCornerElements * rowNum) - (numCornerElements/2);
                if(isOfElement21){
                    a21.push_back(matrixA.at(i));
                    b21.push_back(matrixB.at(i));
                }
                else{
                    a22.push_back(matrixA.at(i));
                    b22.push_back(matrixB.at(i));
                }
            }

            counter++;

            if(counter == numCornerElements){
                counter = 0;
                rowNum++;
            }
        }
    
    /*
        m1 = (a11 + a22) * (b11 + b22);
        m2 = (a21 + a22) * b11;
        m3 = a11 * (b12-b22);
        m4 = a22 * (b21 - b11);
        m5 = (a11 + a12) * b22;
        m6 = (a21 - a11) * (b11 + b12);
        m7 = (a12 - a22) * (b21 + b22);
    */
        m1 = strassenMultiplyRowOrder(addMatrix(a11, a22), addMatrix(b11, b22));
        m2 = strassenMultiplyRowOrder(addMatrix(a21, a22), b11);
        m3 = strassenMultiplyRowOrder(a11, subtractMatrix(b12, b22));
        m4 = strassenMultiplyRowOrder(a22, subtractMatrix(b21, b11));
        m5 = strassenMultiplyRowOrder(addMatrix(a11, a12), b22);
        m6 = strassenMultiplyRowOrder(subtractMatrix(a21,a11), addMatrix(b11, b12));
        m7 = strassenMultiplyRowOrder(subtractMatrix(a12, a22), addMatrix(b21, b22));

    /*
        c11 = m1 + m4 - m5 + m7;
        c12 = m3 + m5;
        c21 = m2 + m4;
        c22 = m1 - m2 + m3 + m6;
    */
        c11 = subtractMatrix(addMatrix(addMatrix(m1, m4), m7), m5);
        c12 = addMatrix(m3, m5);
        c21 = addMatrix(m2, m4);
        c22 = subtractMatrix(addMatrix(addMatrix(m1,m3),m6), m2);
       
       
        //index = (Row x Total Width) + Column
        size_t halfResultWidth = sqrt(numCornerElements);
        size_t resultWidth = halfResultWidth * 2; 

        for (size_t row = 0; row < halfResultWidth; ++row) {
            for (size_t col = 0; col < halfResultWidth; ++col) {
                
                size_t subIdx = row * halfResultWidth + col;

                strassenResult.at(row * resultWidth + col) = c11.at(subIdx);
                strassenResult.at(row * resultWidth + (col + halfResultWidth)) = c12.at(subIdx);
                strassenResult.at((row + halfResultWidth) * resultWidth + col) = c21.at(subIdx);
                strassenResult.at((row + halfResultWidth) * resultWidth + (col + halfResultWidth)) = c22.at(subIdx);
            }
        }
        
    }
    else{

        int a11, a12, a21, a22;
        int b11, b12, b21, b22;
        int m1, m2, m3, m4, m5, m6, m7;
        int c11, c12, c21, c22;

        a11 = matrixA.at(0);
        b11 = matrixB.at(0);
        a12 = matrixA.at(1);
        b12 = matrixB.at(1);
        a21 = matrixA.at(2);
        b21 = matrixB.at(2);
        a22 = matrixA.at(3);
        b22 = matrixB.at(3);

        m1 = (a11 + a22) * (b11 + b22);
        m2 = (a21 + a22) * b11;
        m3 = a11 * (b12-b22);
        m4 = a22 * (b21 - b11);
        m5 = (a11 + a12) * b22;
        m6 = (a21 - a11) * (b11 + b12);
        m7 = (a12 - a22) * (b21 + b22);    
        
        c11 = m1 + m4 - m5 + m7;
        c12 = m3 + m5;
        c21 = m2 + m4;
        c22 = m1 - m2 + m3 + m6;

        strassenResult = {c11, c12, c21, c22};
    }

    return strassenResult;
}