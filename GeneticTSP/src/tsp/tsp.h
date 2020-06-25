//
// Created by francesco on 25/06/20.
//

#ifndef GENETICTSP_SRC_TSP_TSP_H_
#define GENETICTSP_SRC_TSP_TSP_H_

#include "../graph/graph.h"
class TSP {
 public:
  template<typename Key = int>
   void Run(Graph<Key, double> &graph,
           int chromosomeNumber,
           int generationNumber,
           double mutationRate,
           double crossoverRate,
           int workers,
           int seed);
  virtual ~TSP() = default;
};
#endif //GENETICTSP_SRC_TSP_TSP_H_
