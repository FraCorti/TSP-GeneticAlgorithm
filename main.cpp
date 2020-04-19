#include <iostream>
#include <omp.h>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <thread>
#include "gameOfLife.h"




int main(int argc, char *argv[]) {
  int iterationNumber = std::atoi(argv[1]);
  int generatorSeedNumber = std::atoi(argv[2]);
  int rowNumbers = std::atoi(argv[3]);
  int columnNumber = std::atoi(argv[4]);
  int parallelismDegree = 1;
  if(argv[5]){
    parallelismDegree = std::atoi(argv[5]);
  }

  std::cout << "C++ threads " << std::endl;
  for (int threadNumber = 1; threadNumber < 150; threadNumber++) {
    GameOfLife gameOfLife(iterationNumber, generatorSeedNumber, rowNumbers, columnNumber, threadNumber);
    gameOfLife.StandardThreads();
  }

  std::cout << "PRAGMA threads " << std::endl;
  for (int threadNumber = 1; threadNumber < 150; threadNumber++) {
    GameOfLife gameOfLife(iterationNumber, generatorSeedNumber, rowNumbers, columnNumber, threadNumber);
    gameOfLife.PragmaParallel();
  }
  return 0;

}