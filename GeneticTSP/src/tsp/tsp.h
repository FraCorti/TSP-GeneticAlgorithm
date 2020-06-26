//
// Created by francesco on 25/06/20.
//

#ifndef GENETICTSP_SRC_TSP_TSP_H_
#define GENETICTSP_SRC_TSP_TSP_H_

class TSP {
 public:
  virtual ~TSP() = default;
  virtual void Run(int chromosomeNumber,
              int generationNumber,
              double mutationRate,
              double crossoverRate,
              int workers,
              int seed) = 0;
};
#endif //GENETICTSP_SRC_TSP_TSP_H_