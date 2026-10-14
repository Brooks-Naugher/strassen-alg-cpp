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
bool isStrassenMultipliable1D(std::vector<int>& matrix1D){
    
    size_t size = matrix1D.size();
    
    if(size == 0){
        return false;
    }

    //is not a square of 2
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
bool areBothStrassenMultipliable1D(std::vector<int>& matrixA, std::vector<int>& matrixB){
    size_t matrixASize = matrixA.size();
    size_t matrixBSize = matrixB.size();

    bool areSameSize = matrixASize == matrixBSize;
    if(!areSameSize){
        return false;
    }

    bool areStrassenMultipliable = isStrassenMultipliable1D(matrixA) && isStrassenMultipliable1D(matrixB);
    
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
        strassenResult = {matrixA[0]*matrixB[0]};
    }
    else if(numElements != 4){

        std::size_t numCornerElements = numElements / 4;
        
        std::vector<int> a11, a12, a21, a22;
        std::vector<int> b11, b12, b21, b22;
        std::vector<int> m1, m2, m3, m4, m5, m6, m7;
        std::vector<int> c11, c12, c21, c22;

        for(std::size_t i = 0; i < numElements; i++){
            if(i < numCornerElements){
                a11.push_back(matrixA[i]);
                b11.push_back(matrixB[i]);
            } 
            else if(i < numCornerElements * 2){
                a12.push_back(matrixA[i]);
                b12.push_back(matrixB[i]);
            }
            else if(i < numCornerElements * 3){
                a21.push_back(matrixA[i]);
                b21.push_back(matrixB[i]);
            }
            else{
                a22.push_back(matrixA[i]);
                b22.push_back(matrixB[i]); 
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
        m1 = strassenMultiplyCornerOrder(add1dMatrices(a11, a22), add1dMatrices(b11, b22));
        m2 = strassenMultiplyCornerOrder(add1dMatrices(a21, a22), b11);
        m3 = strassenMultiplyCornerOrder(a11, subtract1dMatrices(b12, b22));
        m4 = strassenMultiplyCornerOrder(a22, subtract1dMatrices(b21, b11));
        m5 = strassenMultiplyCornerOrder(add1dMatrices(a11, a12), b22);
        m6 = strassenMultiplyCornerOrder(subtract1dMatrices(a21,a11), add1dMatrices(b11, b12));
        m7 = strassenMultiplyCornerOrder(subtract1dMatrices(a12, a22), add1dMatrices(b21, b22));

    /*
        c11 = m1 + m4 - m5 + m7;
        c12 = m3 + m5;
        c21 = m2 + m4;
        c22 = m1 - m2 + m3 + m6;
    */
        c11 = subtract1dMatrices(add1dMatrices(add1dMatrices(m1, m4), m7), m5);
        c12 = add1dMatrices(m3, m5);
        c21 = add1dMatrices(m2, m4);
        c22 = subtract1dMatrices(add1dMatrices(add1dMatrices(m1,m3),m6), m2);

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

        a11 = matrixA[0];
        b11 = matrixB[0];
        a12 = matrixA[1];
        b12 = matrixB[1];
        a21 = matrixA[2];
        b21 = matrixB[2];
        a22 = matrixA[3];
        b22 = matrixB[3];

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
        strassenResult = {matrixA[0]*matrixB[0]};
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
                    a11.push_back(matrixA[i]);
                    b11.push_back(matrixB[i]);
                } 
                else {
                    a12.push_back(matrixA[i]);
                    b12.push_back(matrixB[i]);                }
            }
            else{
                bool isOfElement21 = i < (numCornerElements * rowNum) - (numCornerElements/2);
                if(isOfElement21){
                    a21.push_back(matrixA[i]);
                    b21.push_back(matrixB[i]);
                }
                else{
                    a22.push_back(matrixA[i]);
                    b22.push_back(matrixB[i]);
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
        m1 = strassenMultiplyRowOrder(add1dMatrices(a11, a22), add1dMatrices(b11, b22));
        m2 = strassenMultiplyRowOrder(add1dMatrices(a21, a22), b11);
        m3 = strassenMultiplyRowOrder(a11, subtract1dMatrices(b12, b22));
        m4 = strassenMultiplyRowOrder(a22, subtract1dMatrices(b21, b11));
        m5 = strassenMultiplyRowOrder(add1dMatrices(a11, a12), b22);
        m6 = strassenMultiplyRowOrder(subtract1dMatrices(a21,a11), add1dMatrices(b11, b12));
        m7 = strassenMultiplyRowOrder(subtract1dMatrices(a12, a22), add1dMatrices(b21, b22));

    /*
        c11 = m1 + m4 - m5 + m7;
        c12 = m3 + m5;
        c21 = m2 + m4;
        c22 = m1 - m2 + m3 + m6;
    */
        c11 = subtract1dMatrices(add1dMatrices(add1dMatrices(m1, m4), m7), m5);
        c12 = add1dMatrices(m3, m5);
        c21 = add1dMatrices(m2, m4);
        c22 = subtract1dMatrices(add1dMatrices(add1dMatrices(m1,m3),m6), m2);
       
       
        //index = (Row x Total Width) + Column
        size_t halfResultWidth = sqrt(numCornerElements);
        size_t resultWidth = halfResultWidth * 2; 

        for (size_t row = 0; row < halfResultWidth; ++row) {
            for (size_t col = 0; col < halfResultWidth; ++col) {
                
                size_t subIdx = row * halfResultWidth + col;

                strassenResult[row * resultWidth + col] = c11[subIdx];
                strassenResult[row * resultWidth + (col + halfResultWidth)] = c12[subIdx];
                strassenResult[(row + halfResultWidth) * resultWidth + col] = c21[subIdx];
                strassenResult[(row + halfResultWidth) * resultWidth + (col + halfResultWidth)] = c22[subIdx];
            }
        }
        
    }
    else{

        int a11, a12, a21, a22;
        int b11, b12, b21, b22;
        int m1, m2, m3, m4, m5, m6, m7;
        int c11, c12, c21, c22;

        a11 = matrixA[0];
        b11 = matrixB[0];
        a12 = matrixA[1];
        b12 = matrixB[1];
        a21 = matrixA[2];
        b21 = matrixB[2];
        a22 = matrixA[3];
        b22 = matrixB[3];

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
 * @brief Multiplies two matrices using Strassen Algorithm
 * 
 * Multiplies two same-sized, Strassen multipliable
 * matrices (2d Array). The matrices are in "Row Order".
 * Row Order Matrix  looks like this: If
        m = [ [1,2,3,4], [5,6,7,8], [9,10,11,12], [13,14,15,16] ]
        Each Row: 
        m1 = 1,2,3,4 | 
        m2 = 4,6,7,8 | 
        m3 = 9,10,11,12 | 
        m4 = 13,14,15,16
    @param matrixA
    @param matrixB
    @return matrixC - multiplied result of matrixA and matrixB

 */
std::vector<std::vector<int>> strassenMultiply2DVector(const std::vector<std::vector<int>>& matrixA, const std::vector<std::vector<int>>& matrixB){

    size_t numRows = matrixA.size();
    size_t numElements = numRows * numRows;
    std::vector<std::vector<int>> strassenResult;

    if(numElements == 1){
        strassenResult = {{matrixA[0][0]*matrixB[0][0]}};
    }
    else if(numElements != 4){

        std::vector<std::vector<int>> a11, a12, a21, a22;
        std::vector<std::vector<int>> b11, b12, b21, b22;
        std::vector<std::vector<int>> m1, m2, m3, m4, m5, m6, m7;
        std::vector<std::vector<int>> c11, c12, c21, c22;

        std::size_t rowSize = matrixA[0].size();


        for(size_t i =  0; i < numRows; i++){

            std::vector<int> aFirst, aSecond, bFirst, bSecond;

            for(std::size_t j = 0; j < rowSize; j++){

                if(j < rowSize / 2){
                    aFirst.push_back(matrixA[i][j]);
                    bFirst.push_back(matrixB[i][j]);
                }
                else{
                    aSecond.push_back(matrixA[i][j]);
                    bSecond.push_back(matrixB[i][j]);
                }
            }


            if(i < numRows / 2){
                a11.push_back(aFirst);
                a12.push_back(aSecond);
                b11.push_back(bFirst);
                b12.push_back(bSecond);
            } 
            else{
                a21.push_back(aFirst);
                a22.push_back(aSecond);
                b21.push_back(bFirst);
                b22.push_back(bSecond);            
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
            m1 = strassenMultiply2DVector(add2dMatrices(a11, a22), add2dMatrices(b11, b22));
            m2 = strassenMultiply2DVector(add2dMatrices(a21, a22), b11);
            m3 = strassenMultiply2DVector(a11, subtract2dMatrices(b12, b22));
            m4 = strassenMultiply2DVector(a22, subtract2dMatrices(b21, b11));
            m5 = strassenMultiply2DVector(add2dMatrices(a11, a12), b22);
            m6 = strassenMultiply2DVector(subtract2dMatrices(a21, a11), add2dMatrices(b11, b12));
            m7 = strassenMultiply2DVector(subtract2dMatrices(a12, a22), add2dMatrices(b21, b22));


            /*
                c11 = m1 + m4 - m5 + m7;
                c12 = m3 + m5;
                c21 = m2 + m4;
                c22 = m1 - m2 + m3 + m6;
            */

            c11 = subtract2dMatrices(add2dMatrices(m1, add2dMatrices(m4, m7)), m5);
            c12 = add2dMatrices(m3, m5);
            c21 = add2dMatrices(m2, m4);
            c22 = add2dMatrices(subtract2dMatrices(m1, m2), add2dMatrices(m3, m6));
            
            std::size_t cornerSideSize = c11.size();

            for(size_t i = 0; i < cornerSideSize; i++){
                std::vector<int> row(rowSize, 0);
                for(size_t j = 0; j < cornerSideSize; j++){
                    row[j] = c11[i][j];
                    row[j + cornerSideSize] = c12[i][j];
                }
                strassenResult.push_back(row);
            }

            for(size_t i = 0; i < cornerSideSize; i++){
                std::vector<int> row(rowSize, 0);
                for(size_t j = 0; j < cornerSideSize; j++){
                    row[j] = c21[i][j];
                    row[j + cornerSideSize] = c22[i][j];
                }
                strassenResult.push_back(row);
            }

    }
    else{

        int a11, a12, a21, a22;
        int b11, b12, b21, b22;
        int m1, m2, m3, m4, m5, m6, m7;
        int c11, c12, c21, c22;

        a11 = matrixA[0][0];
        b11 = matrixB[0][0];
        a12 = matrixA[0][1];
        b12 = matrixB[0][1];
        a21 = matrixA[1][0];
        b21 = matrixB[1][0];
        a22 = matrixA[1][1];
        b22 = matrixB[1][1];

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

        strassenResult = {{c11, c12}, {c21, c22}};

    }

    return strassenResult;
}