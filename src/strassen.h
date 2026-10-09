#include <vector>
#include <cmath>
#include "matrix_operations.h"

//https://en.wikipedia.org/wiki/Strassen_algorithm

/**
 * @brief Multiplies two matrices using Strassen Algorithm.
 * 
 * Multiplies two same-sized, Strassen multipliable
 * matrices. The matrices are in "Corner Order".
 * Corner Order Matrix  looks like this: If
        m = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
        Then: 
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

    int numElements = matrixA.size();
    std::vector<int> strassenResult;

    if(numElements == 1){
        strassenResult = {matrixA.at(0)*matrixB.at(0)};
    }
    else if(numElements != 4){

        int numCornerElements = numElements / 4;
        
        std::vector<int> a11, a12, a21, a22;
        std::vector<int> b11, b12, b21, b22;
        std::vector<int> m1, m2, m3, m4, m5, m6, m7;
        std::vector<int> c11, c12, c21, c22;

        for(int i = 0; i < numElements; i++){
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

std::vector<int> strassenMultiplyRowOrder(const std::vector<int>& matrixA, const std::vector<int>& matrixB){
    

    int numElements = matrixA.size();
    std::vector<int> strassenResult;

    if(numElements == 1){
        strassenResult = {matrixA.at(0)*matrixB.at(0)};
    }
    else if(numElements != 4){

        int numCornerElements = numElements / 4;
        int numRowElementsNeeded = numCornerElements / 2;
        std::vector<int> a11, a12, a21, a22;
        std::vector<int> b11, b12, b21, b22;
        std::vector<int> m1, m2, m3, m4, m5, m6, m7;
        std::vector<int> c11, c12, c21, c22;
        
        /*
            v = {1,2,3,4  5,6,7,8  9,10,11,12  13,14,15,16}
                 0 1      4 5      8           12
                
            v = {1,2,3,4,5,6,7,8  1,2,3,4,5,6,7,8  1,2,3,4,5,6,7,8  1,2,3,4,5,6,7,8  1,2,3,4,5,6,7,8  1,2,3,4,5,6,7,8  1,2,3,4,5,6,7,8  1,2,3,4,5,6,7,8}
        */     //a11 = 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4

        
    /*
        for(int i = 0; i < numElements; i++){

            if(i < numElements/2){
                bool isWithFirstRowElement11 = i < numCornerElements && i < numCornerElements / 2;
                bool isWithSecondRowElemenet11 = i > numCornerElements && i < (numCornerElements*2 - numRowElementsNeeded );

                bool 
                if(isWithFirstRowElement11 || isWithSecondRowElemenet11){
                    a11.push_back(matrixA.at(i));
                    b11.push_back(matrixB.at(i));
                } else{
                    a12.push_back(matrixA.at(i));
                    b12.push_back(matrixB.at(i));
                }
            } else{
                bool isWithFirstRowA21 
                bool isWithSecondRowA22
                if(){

                } else{
                    
                }
            }
        }
    */

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