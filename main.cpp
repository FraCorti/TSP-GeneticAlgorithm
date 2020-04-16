#include <iostream>
#include <omp.h>
#include <vector>
#include <cstdlib>
#include <chrono>

int main(int argc, char *argv[]) {
  int iterationNumber = std::atoi(argv[1]);
  int generatorSeedNumber = std::atoi(argv[2]);
  int rowNumbers = std::atoi(argv[3]);
  int columnNumber = std::atoi(argv[4]);
  int parallelismDegree = 1;
  if(argv[5]){
    parallelismDegree = std::atoi(argv[5]);
  }

  srand(generatorSeedNumber);

  //! initialize the two matrix to zero
  int* gameOfLifeStart = new int[rowNumbers*columnNumber]{};
  int* gameOfLifeEnd = new int[rowNumbers*columnNumber]{};

  //! fill up start matrix
  for(int row = 1; row < rowNumbers -1; row++){
    for(int column = 1; column < columnNumber - 1; column++){
      gameOfLifeStart[row*columnNumber +column] = rand() % 2;
    }
  }

  auto start = std::chrono::system_clock::now();

  for(int iteration = 0; iteration < iterationNumber; iteration++) {

    //! compute the number of alive neighbours
    for(int row = 1; row < rowNumbers-1; row++) {

      //! sequential code + flag for openmp parallel for
      // #pragma omp parallel for num_threads(parallelismDegree)
      for (int column = 1; column < columnNumber - 1; column++) {
        int counter = gameOfLifeStart[(row*columnNumber+column) -(columnNumber+1)] +
            gameOfLifeStart[(row*columnNumber+column) -columnNumber]  +
            gameOfLifeStart[(row*columnNumber+column) -(columnNumber-1)]  +
            gameOfLifeStart[(row*columnNumber+column) -1]  +
            gameOfLifeStart[(row*columnNumber+column) +1]  +
            gameOfLifeStart[(row*columnNumber+column) +(columnNumber-1)]  +
            gameOfLifeStart[(row*columnNumber+column) +columnNumber]  +
            gameOfLifeStart[(row*columnNumber+column) +(columnNumber+1)];
        if(counter >= 2 && counter <= 3){
          // an empty cell with 2 neighbours shouldn't spawn an individual
          gameOfLifeEnd[row*columnNumber +column] = 1;
        } else{
          gameOfLifeEnd[row*columnNumber +column] = 0;
        }
      }
    }
    //! print current configuration
    for (int row = 0; row < rowNumbers ; row++) {
      for (int column = 0; column < columnNumber; column++) {
        if(gameOfLifeStart[row*columnNumber+column] == 1){
          std::cout << "O";
        } else{
          std::cout<< "_";
        }
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    //! swap the arrays
    std::swap(gameOfLifeStart, gameOfLifeEnd);
  }

  auto end = std::chrono::system_clock::now();
  std::cout<< "Time: "<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()<<"ms" << std::endl;
  return 0;
}