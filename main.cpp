#include <iostream>
#include <omp.h>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <thread>
#include "gameOfLife.h"


int main(int argc, char *argv[]) {
  if (argc == 1) {
    std::cout << "Usage is " << argv[0]
              << " numberOfIteration seed N M nw" << std::endl;
    return (-1);
  }

  int iterationNumber = std::atoi(argv[1]);
  int generatorSeedNumber = std::atoi(argv[2]);
  int rowNumbers = std::atoi(argv[3]);
  int columnNumber = std::atoi(argv[4]);
  int parallelismDegree = 1;
  if(argv[5]){
    parallelismDegree = std::atoi(argv[5]);
  }

  GameOfLife gameOfLife(iterationNumber, generatorSeedNumber, rowNumbers, columnNumber, parallelismDegree);

  //! Sequential
  gameOfLife.Sequential();

  //! Standard
  gameOfLife.StandardThreads();

  //! OpenMP
  gameOfLife.OpenMP();

  return 0;
}