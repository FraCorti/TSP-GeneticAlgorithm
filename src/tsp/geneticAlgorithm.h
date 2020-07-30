//
// Created by francesco on 25/06/20.
//

#ifndef GENETICTSP_SRC_TSP_GENETICALGORITHM_H_
#define GENETICTSP_SRC_TSP_GENETICALGORITHM_H_

#include <iostream>
#include <vector>
#include <functional>
#include <chrono>

//! let the user decide the precision of intermediate results
using precision = double;

class GeneticAlgorithm {
 public:
  virtual ~GeneticAlgorithm() = default;
  virtual void Run(int chromosomeNumber,
                   int generationNumber,
                   double mutationRate,
                   double crossoverRate,
                   int workers,
                   int seed) = 0;
};
#endif //GENETICTSP_SRC_TSP_GENETICALGORITHM_H_