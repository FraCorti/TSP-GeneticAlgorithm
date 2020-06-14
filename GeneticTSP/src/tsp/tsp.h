//
// Created by francesco on 31/05/20.
//

#ifndef GENETICTSP_SRC_TSP_TSP_H_
#define GENETICTSP_SRC_TSP_TSP_H_

#include <vector>
#include "../graph/graph.h"
#include <functional>

template<typename Key = int>
class TSP {
 private:
  //! vector containing pairs of <fitnessScore, chromosome>
  std::vector<std::pair<double, std::vector<std::pair<Key, double>>>> population;
  std::random_device rd;    // Will be used to obtain a seed for the random number engine
  std::mt19937 gen{rd()};   //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<double> unif{0, 1};

  const int generationNumber;
  const int workers;
  // TODO: is inline useful in template method?
  inline double generateProbability();
  inline void evaluate(Graph<Key, double> &graph);
  inline void crossover();
  inline void fitness();
  inline void printPopulation();
  inline void selection();
  inline void sortProbabilities();
  inline void generate();
  static bool CompareProbability(const std::pair<Key, double> &a, const std::pair<Key, double> &b);
  static bool CompareChromosome(const std::pair<double, std::vector<std::pair<Key, double>>> &a,
                                const std::pair<double, std::vector<std::pair<Key, double>>> &b);
 public:
  explicit TSP(Graph<Key, double> &graph, int generationNumber = 500, int workers = 1);
};

template<typename Key>
TSP<Key>::TSP(Graph<Key, double> &graph, int generationNumber, int workers) :
    generationNumber(generationNumber), workers(workers), population(graph.GetNodesNumber()) {

  int nodesNumber = graph.GetNodesNumber();

  //! generate first population
  std::for_each(population.begin(),
                population.end(),
                [&](std::pair<double, std::vector<std::pair<Key, double>>> &chromosome) {
                  chromosome.second.resize(nodesNumber);
                  auto i = graph.GetMapIterator();

                  //! create chromosome probabilities
                  std::for_each(chromosome.second.begin(),
                                chromosome.second.end(),
                                [&](std::pair<Key, double> &cell) {
                                  cell.first = i->first;
                                  cell.second = generateProbability();
                                  i++;
                                });
                });

  //! start genetic algorithm
  for (int generation = 1; generation <= generationNumber; generation++) {
    sortProbabilities();
    evaluate(graph);
    fitness();
    selection();
    crossover();
    generate();
  }
  printPopulation();
}

/*** For each chromosome compute the fitness score
 *   by summing all the edges value then sort population vector
 *   over the fitness score value
 */
template<typename Key>
void TSP<Key>::evaluate(Graph<Key, double> &graph) {

  int chromosomeSize = population.begin()->second.size();
  std::for_each(population.begin(),
                population.end(),
                [&](std::pair<double, std::vector<std::pair<Key, double>>> &chromosome) {
                  double chromosomeScore = 0;

                  // TODO: can be done with two iterator (maybe is better?)
                  for (int i = 0; i < chromosomeSize - 1; i++) {
                    double edgeValue = graph.GetEdgeValue(chromosome.second[i].first, chromosome.second[i + 1].first);
                    chromosomeScore += edgeValue;
                  }

                  //! retrieve last edge (first and last chromosome element) value and store fitness score
                  chromosome.first += (chromosomeScore
                      + graph.GetEdgeValue(chromosome.second[0].first, chromosome.second[chromosomeSize - 1].first));
                });
  //std::sort(population.begin(), population.end(), CompareChromosome);
}

/***
 *
 */
template<typename Key>
void TSP<Key>::generate() {

}

/***
 *
 */
template<typename Key>
void TSP<Key>::crossover() {

}

/*** Select the best
 *
 */
template<typename Key>
void TSP<Key>::selection() {

}
template<typename Key>
double TSP<Key>::generateProbability() {
  return unif(gen);
}

/*** Iterate over the chromosomes and sorting in descending order
 *   over to the probability
 */
template<typename Key>
void TSP<Key>::sortProbabilities() {
  std::for_each(population.begin(),
                population.end(),
                [&](std::pair<double, std::vector<std::pair<Key, double>>> &chromosome) {
                  std::sort(chromosome.second.begin(), chromosome.second.end(), CompareProbability);
                });
}

template<typename Key>
void TSP<Key>::printPopulation() {
  //! print current population
  std::for_each(population.begin(),
                population.end(),
                [&](std::pair<double, std::vector<std::pair<Key, double>>> &chromosome) {
                  std::cout << "Fitness Score: " << chromosome.first << std::endl;
                  std::for_each(chromosome.second.begin(),
                                chromosome.second.end(),
                                [&](std::pair<Key, double> &cell) {
                                  std::cout << cell.first << " " << cell.second << std::endl;
                                });
                  std::cout << std::endl;
                });
}

/*** Compare two pair of the form <nodeKey, propability>
 *
 * @return
 */
template<typename Key>
bool TSP<Key>::CompareProbability(const std::pair<Key, double> &a, const std::pair<Key, double> &b) {
  return (a.second > b.second);
}

template<typename Key>
bool TSP<Key>::CompareChromosome(const std::pair<double, std::vector<std::pair<Key, double>>> &a,
                                 const std::pair<double, std::vector<std::pair<Key, double>>> &b) {
  return (a.first > b.first);
}
template<typename Key>
void TSP<Key>::fitness() {

}
#endif //GENETICTSP_SRC_TSP_TSP_H_