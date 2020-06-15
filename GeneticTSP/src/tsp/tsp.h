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

  //! vector of chromosome
  std::vector<std::vector<std::pair<Key, double>>> population;

  std::random_device rd;    // Will be used to obtain a seed for the random number engine
  std::mt19937 gen{rd()};   //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<double> unif{0, 1};

  const int generationNumber;
  const int workers;
  // TODO: is inline useful in template method?
  inline double generateProbability();
  inline void evaluate(Graph<Key, double> &graph,
                       std::vector<std::pair<double, int>> &chromosomesScoreIndex,
                       double &evaluationsAverage);
  inline void crossover();
  inline void fitness(std::vector<std::pair<double, int>> &chromosomesScoreIndex, double evaluationsAverage);
  inline void printPopulation();
  inline void selection(std::vector<std::pair<double, int>> &chromosomesProbabilityIndex,
                        std::vector<std::vector<std::pair<Key, double>>> &intermediatePopulation);
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
                [&](std::vector<std::pair<Key, double>> &chromosome) {
                  chromosome.resize(nodesNumber);
                  auto i = graph.GetMapIterator();

                  //! create chromosome probabilities
                  std::for_each(chromosome.begin(),
                                chromosome.end(),
                                [&](std::pair<Key, double> &cell) {
                                  cell.first = i->first;
                                  cell.second = generateProbability();
                                  i++;
                                });
                });

  //! pairs of <fitness score, chromosomeIndex>
  std::vector<std::pair<double, int>> chromosomesScoreIndex(nodesNumber);

  std::vector<std::vector<std::pair<Key, double>>> intermediatePopulation(population.size());
  double evaluationsAverage = 0;

  //! start genetic algorithm
  for (int generation = 1; generation <= generationNumber; generation++) {
    sortProbabilities();
    evaluate(graph, chromosomesScoreIndex, evaluationsAverage);
    fitness(chromosomesScoreIndex, evaluationsAverage);
    selection(chromosomesScoreIndex, intermediatePopulation);
    crossover();
    generate();

    evaluationsAverage = 0;
    std::swap(intermediatePopulation, population);  //TODO: check if works properly
  }
  printPopulation();
}

/*** For each chromosome compute the fitness score
 *   by summing all the edges value then sort population vector
 *   over the fitness score value
 */
template<typename Key>
void TSP<Key>::evaluate(Graph<Key, double> &graph,
                        std::vector<std::pair<double, int>> &chromosomesScoreIndex,
                        double &evaluationsAverage) {

  int chromosomeSize = population.begin()->size();
  int iteratorIndex = 0;
  std::for_each(population.begin(),
                population.end(),
                [&](std::vector<std::pair<Key, double>> &chromosome) {
                  double chromosomeScore = 0;

                  // TODO: can be done with two iterator (better?)
                  for (int i = 0; i < chromosomeSize - 1; i++) {
                    double edgeValue = graph.GetEdgeValue(chromosome[i].first, chromosome[i + 1].first);
                    chromosomeScore += edgeValue;
                  }

                  //! retrieve last edge value and store fitness score
                  chromosomesScoreIndex[iteratorIndex].first += (chromosomeScore
                      + graph.GetEdgeValue(chromosome[0].first, chromosome[chromosomeSize - 1].first));
                  chromosomesScoreIndex[iteratorIndex].second = iteratorIndex;

                  evaluationsAverage += chromosomesScoreIndex[iteratorIndex].first;
                  iteratorIndex++;
                });

  evaluationsAverage /= chromosomesScoreIndex.size();
}

/*** Compute reproductive opportunity for each pair of chromosomesScoreIndex.
 *   This is done dividing the score (in our case the length of the path)
 *   by the average score obtained
 *
 */
template<typename Key>
void TSP<Key>::fitness(std::vector<std::pair<double, int>> &chromosomesScoreIndex, double evaluationsAverage) {
  /*
  //! print chromosome index with fitness value for test purpose
  std::for_each(chromosomesScoreIndex.begin(),
                chromosomesScoreIndex.end(), [&](std::pair<double, int> &currentIndex) {
        std::cout << "Evaluation: " << currentIndex.first << " Index: " << currentIndex.second << std::endl;
      }); */

  //! for each chromosome compute and store probability of survive
  std::for_each(chromosomesScoreIndex.begin(),
                chromosomesScoreIndex.end(),
                [&](std::pair<double, int> &chromosome) {
                  //! formula adapted from "A genetic algorithm tutorial" Darrel Whitley
                  chromosome.first = evaluationsAverage / chromosome.first;
                });

  //! sort obtained probability in increasing order (ordine crescente)
  std::sort(chromosomesScoreIndex.begin(),
            chromosomesScoreIndex.end(),
            [&](const std::pair<double, int> &a, const std::pair<double, int> &b) {
              return a.first < b.first;
            });

  /*
  //! print chromosome index with probability for test purpose
  std::for_each(chromosomesScoreIndex.begin(),
                chromosomesScoreIndex.end(), [&](std::pair<double, int> &currentIndex) {
        std::cout << "Probability: " << currentIndex.first << " Index: " << currentIndex.second << std::endl;
      }); */
}

/*** Produce the intermediate population based on the probability of reproduction computed before then
 *   the new intermediate population is saved inside population
 *
 *   @param chromosomesProbabilityIndex Vector of <reproduceProbability, chromosomeIndex> sorted in increasing order
 */
template<typename Key>
void TSP<Key>::selection(std::vector<std::pair<double, int>> &chromosomesProbabilityIndex,
                         std::vector<std::vector<std::pair<Key, double>>> &intermediatePopulation) {

  double probabilitiesSum = 0;
  std::vector<double> randomProbabilities(chromosomesProbabilityIndex.size());
  std::vector<int> indexNewPopulation;
  indexNewPopulation.reserve(chromosomesProbabilityIndex.size());

  //! divide the probability space between the chromosome over their reproduction factor
  std::for_each(chromosomesProbabilityIndex.begin(),
                chromosomesProbabilityIndex.end(),
                [&](std::pair<double, int> &chromosome) {
                  probabilitiesSum += chromosome.first;
                  chromosome.first = probabilitiesSum;
                });
  /*
  //! print chromosome index with probability for test purpose
  std::for_each(chromosomesProbabilityIndex.begin(),
                chromosomesProbabilityIndex.end(), [&](std::pair<double, int> &currentIndex) {
        std::cout << "Probability: " << currentIndex.first << " Index: " << currentIndex.second << std::endl;
      }); */

  //! generate random numbers inside probability space and sort them
  std::uniform_real_distribution<double> probabilitySpaceInterval(0, probabilitiesSum);
  std::for_each(randomProbabilities.begin(),
                randomProbabilities.end(), [&](double &currentProbability) {
        currentProbability = probabilitySpaceInterval(gen);
      });
  std::sort(randomProbabilities.begin(), randomProbabilities.end());

  //! fill the vector of new population indexes
  auto randomProbabilitiesIt = randomProbabilities.begin();
  for (std::pair<double, int> &chromosome: chromosomesProbabilityIndex) {
    for (; randomProbabilitiesIt != randomProbabilities.end() && *randomProbabilitiesIt <= chromosome.first;
           randomProbabilitiesIt++) {
      indexNewPopulation.push_back(chromosome.second);
    }
  }

  /*
 //! print new population index
 for(int & a : indexNewPopulation){
   std::cout << "Index "<< a << std::endl;
 }

 std::cout << "population before swap "<< std::endl;
 printPopulation(); */

  //! create intermediate population
  auto populationIter = intermediatePopulation.begin();
  std::for_each(indexNewPopulation.begin(),
                indexNewPopulation.end(), [&](int &index) {
        *populationIter = population[index];
        populationIter++;
      });

  /*
  std::cout << "population after swap "<< std::endl;
  printPopulation(); */

  /*
  //! print generated probability for test purpose
  std::for_each(randomProbabilities.begin(),
                randomProbabilities.end(), [&](double &currentProbability) {
        std::cout << currentProbability << std::endl;
      }); */

  /*
  //! print chromosome index with increasing probability for test purpose
  std::for_each(chromosomesProbabilityIndex.begin(),
                chromosomesProbabilityIndex.end(), [&](std::pair<double, int> &currentIndex) {
        std::cout << "Probability: " << currentIndex.first << " Index: " << currentIndex.second << std::endl;
      }); */
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
                [&](std::vector<std::pair<Key, double>> &chromosome) {
                  std::sort(chromosome.begin(),
                            chromosome.end(),
                            [&](const std::pair<Key, double> &a, const std::pair<Key, double> &b) {
                              return a.second > b.second;
                            });
                });
}

template<typename Key>
void TSP<Key>::printPopulation() {
  //! print current population
  std::for_each(population.begin(),
                population.end(),
                [&](std::vector<std::pair<Key, double>> &chromosome) {
                  std::for_each(chromosome.begin(),
                                chromosome.end(),
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
#endif //GENETICTSP_SRC_TSP_TSP_H_