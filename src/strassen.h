#include <vector>

//https://en.wikipedia.org/wiki/Strassen_algorithm

//Divide and Conquer

std::vector<std::vector<short>> strassenMultiply(std::vector<std::vector<short>> matrixA, std::vector<std::vector<short>> matrixB){

    short numCornerElements = (matrixA.size() / 2) * (matrixB.size() / 2);

    //matrix = { {1, 2}, {3, 4} }
    if(numCornerElements == 1){
        short a11, a12, a21, a22, b11, b12, b21, b22; 
        short m1, m2, m3, m4, m5, m6, m7;
        short c11, c12, c21, c22;
        std::vector<std::vector<short>> strassenResult;

        a11 = matrixA.at(0).at(0);
        a12 = matrixA.at(0).at(1);
        a21 = matrixA.at(1).at(0);
        a22 = matrixA.at(1).at(1);

        b11 = matrixB.at(0).at(0);
        b12 = matrixB.at(0).at(1);
        b21 = matrixB.at(1).at(0);
        b22 = matrixB.at(1).at(1);

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

        strassenResult = {{c11, c12},{c21, c22}};
        return strassenResult;
    }

}