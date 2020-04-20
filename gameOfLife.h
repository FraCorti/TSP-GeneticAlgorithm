//
// Created by checco on 17/04/20.
//

#ifndef FIRSTASSIGNMENTPARALLEL__GAMEOFLIFE_H_
#define FIRSTASSIGNMENTPARALLEL__GAMEOFLIFE_H_

#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>
#include <cmath>
#include <deque>

class GameOfLife {
 private:
  const int iterationNumber;
  const int generatorSeedNumber;
  const int rowNumbers;
  const int columnsNumbers;
  const int parallelismDegree;
  bool *gameOfLifeStart;
  bool *gameOfLifeEnd;
  inline void PrintCurrentEpoch();
  inline void Swap();
  static void WorkerComputation(bool *gameOfLifeStart,
                                bool *gameOfLifeEnd,
                                const int indexStart,
                                const int indexEnd,
                                const int columns);
 public:
  GameOfLife(int iteration_number,
             int generator_SeedNumber,
             int row_numbers,
             int columns_numbers,
             int parallelism_degree);
  void Sequential();
  void StandardThreads();
  void OpenMP();
};

GameOfLife::GameOfLife(int iteration_number,
                       int generator_SeedNumber,
                       int row_numbers,
                       int columns_numbers,
                       int parallelism_degree
) : iterationNumber(iteration_number),
    generatorSeedNumber(generator_SeedNumber),
    rowNumbers(row_numbers),
    columnsNumbers(columns_numbers),
    parallelismDegree(parallelism_degree) {
  std::srand(generatorSeedNumber);
  gameOfLifeStart = new bool[rowNumbers * columnsNumbers]{};
  gameOfLifeEnd = new bool[rowNumbers * columnsNumbers]{};

  //! fill up start matrix
  for (int row = 1; row < rowNumbers - 1; row++) {
    for (int column = 1; column < columnsNumbers - 1; column++) {
      gameOfLifeStart[row * columnsNumbers + column] = rand() % 2;
    }
  }
  PrintCurrentEpoch();
}

//! sequential version
void GameOfLife::Sequential() {
  auto start = std::chrono::system_clock::now();
  for (int iteration = 0; iteration < iterationNumber; iteration++) {
    for (int row = 1; row < columnsNumbers - 1; row++) {
      for (int column = 1; column < rowNumbers - 1; column++) {
        const int counter = gameOfLifeStart[(row - 1) * columnsNumbers + column - 1] +
            gameOfLifeStart[(row - 1) * columnsNumbers + column] +
            gameOfLifeStart[(row - 1) * columnsNumbers + column + 1] +
            gameOfLifeStart[(row) * columnsNumbers + column - 1] +
            gameOfLifeStart[(row) * columnsNumbers + column + 1] +
            gameOfLifeStart[(row + 1) * columnsNumbers + column - 1] +
            gameOfLifeStart[(row + 1) * columnsNumbers + column] +
            gameOfLifeStart[(row + 1) * columnsNumbers + column + 1];
        if (counter == 3) {
          gameOfLifeEnd[(row * columnsNumbers + column)] = true;
        } else
          gameOfLifeEnd[(row * columnsNumbers + column)] =
              counter == 2 && gameOfLifeStart[(row * columnsNumbers + column)];
      }
    }
    Swap();
    PrintCurrentEpoch();
  }
  auto end = std::chrono::system_clock::now();
  std::cout << "Sequential time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << "ms"
            << std::endl;
}

//! single C++ thread computation
void GameOfLife::WorkerComputation(bool *gameOfLifeStart,
                                   bool *gameOfLifeEnd,
                                   const int indexStart,
                                   const int indexEnd,
                                   const int columns) {

  for (int currentIndex = indexStart; currentIndex <= indexEnd; currentIndex++) {

    // exclude matrix border
    if (!(currentIndex % columns) || !((currentIndex + 1) % columns))
      continue;

    const int counter = gameOfLifeStart[currentIndex - (columns + 1)] +
        gameOfLifeStart[currentIndex - (columns)] +
        gameOfLifeStart[currentIndex - (columns - 1)] +
        gameOfLifeStart[currentIndex - 1] +
        gameOfLifeStart[currentIndex + 1] +
        gameOfLifeStart[currentIndex + (columns - 1)] +
        gameOfLifeStart[currentIndex + columns] +
        gameOfLifeStart[currentIndex + (columns + 1)];
    if (counter == 3) {
      gameOfLifeEnd[currentIndex] = true;
    } else
      gameOfLifeEnd[currentIndex] =
          counter == 2 && gameOfLifeStart[currentIndex];
  }
}

//! parallel version with C++ threads
void GameOfLife::StandardThreads() {
  auto start = std::chrono::system_clock::now();
  std::vector<std::thread> computationThreads(parallelismDegree);
  int workerCellsToCompute = std::floor(static_cast<float >(rowNumbers - 2) * static_cast<float>(columnsNumbers)
                                            / static_cast<float>(parallelismDegree));

  for (int iteration = 0; iteration < iterationNumber; iteration++) {
    int remainedElements = ((rowNumbers - 2) * columnsNumbers) % parallelismDegree;
    auto computationThreadIterator = computationThreads.begin();
    int indexStart = columnsNumbers;
    for (int i = 0; i < parallelismDegree; i++) {
      int indexEnd = indexStart + (workerCellsToCompute - 1);

      // consider extra elements
      if (remainedElements) {
        indexEnd++;
        remainedElements--;
      }

      computationThreads.emplace_back(std::thread(
          &GameOfLife::WorkerComputation,
          gameOfLifeStart,
          gameOfLifeEnd,
          indexStart,
          indexEnd,
          columnsNumbers)
      );
      indexStart = indexEnd + 1;
      computationThreadIterator++;
    }

    for (std::thread &thread :computationThreads) {
      if (thread.joinable())
        thread.join();
    }
    computationThreads.clear();
    Swap();
    PrintCurrentEpoch();
  }
  auto end = std::chrono::system_clock::now();
  std::cout << parallelismDegree << " " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << std::endl;

}

//! parallel version using OpenMP
void GameOfLife::OpenMP() {

  auto start = std::chrono::system_clock::now();
  for (int iteration = 0; iteration < iterationNumber; iteration++) {
#pragma omp parallel for num_threads(parallelismDegree)
    for (int row = 1; row < rowNumbers - 1; row++) {
#pragma omp parallel for num_threads(parallelismDegree)
      for (int column = 1; column < columnsNumbers - 1; column++) {
        int counter = gameOfLifeStart[(row - 1) * columnsNumbers + column - 1] +
            gameOfLifeStart[(row - 1) * columnsNumbers + column] +
            gameOfLifeStart[(row - 1) * columnsNumbers + column + 1] +
            gameOfLifeStart[(row) * columnsNumbers + column - 1] +
            gameOfLifeStart[(row) * columnsNumbers + column + 1] +
            gameOfLifeStart[(row + 1) * columnsNumbers + column - 1] +
            gameOfLifeStart[(row + 1) * columnsNumbers + column] +
            gameOfLifeStart[(row + 1) * columnsNumbers + column + 1];
        if (counter == 3) {
          gameOfLifeEnd[(row * columnsNumbers + column)] = true;
        } else
          gameOfLifeEnd[(row * columnsNumbers + column)] =
              counter == 2 && gameOfLifeStart[(row * columnsNumbers + column)];
      }
    }
    Swap();
    PrintCurrentEpoch();
  }
  auto end = std::chrono::system_clock::now();
  std::cout << parallelismDegree << " " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << std::endl;
}

//! print current configuration of the board
void GameOfLife::PrintCurrentEpoch() {
  for (int row = 0; row < rowNumbers; row++) {
    for (int column = 0; column < columnsNumbers; column++) {
      std::cout << gameOfLifeStart[row * columnsNumbers + column] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

//! swap the arrays
void GameOfLife::Swap() {
  std::swap(gameOfLifeEnd, gameOfLifeStart);
}

#endif //FIRSTASSIGNMENTPARALLEL__GAMEOFLIFE_H_