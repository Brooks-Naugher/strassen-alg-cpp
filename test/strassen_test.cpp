#include <iostream>
#include <vector>
#include "../src/strassen.h"

void strassenTest(){
    std::vector<std::vector<short>> matrixA = {{1,2},{3,4}};
    std::vector<std::vector<short>> matrixB = {{5,6},{7,8}};

    std::vector<std::vector<short>> matrixC = strassenMultiply(matrixA, matrixB);

    for(std::vector<short> row : matrixC){
        std::cout<<"[ ";
        for(short element : row){
            std::cout<<element<<" ";
        }
       std::cout<<"]\n";
    }

}

int main(){
    strassenTest();
    return 0;
}