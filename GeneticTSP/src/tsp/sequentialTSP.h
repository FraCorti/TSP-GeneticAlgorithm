//
// Created by francesco on 31/05/20.
//

#ifndef GENETICTSP_SRC_TSP_SEQUENTIALTSP_H_
#define GENETICTSP_SRC_TSP_SEQUENTIALTSP_H_

#include <vector>
#include "tsp.h"
#include <functional>
#include <chrono>

template<typename Key = int>
class TSPSequential : public TSP{ //
 private:

  //! vector of chromosomes
  std::vector<std::vector<Key>> population;

  std::random_device rd;    // Will be used to obtain a seed for the random number engine
  std::mt19937 gen{rd()};   //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<double> unif{0, 1};
  Graph<Key, double>& graph;
  inline void evaluate(Graph<Key, double> &graph,
                       std::vector<std::pair<double, int>> &chromosomesScoreIndex,
                       double &evaluationsAverage);
  inline void mutation(std::vector<std::vector<Key>> &intermediatePopulation, double mutationRate);
  inline void crossover(std::vector<std::vector<Key>> &intermediatePopulation, double crossoverRate);
  inline void fitness(std::vector<std::pair<double, int>> &chromosomesScoreIndex, double evaluationsAverage);
  inline void selection(std::vector<std::pair<double, int>> &chromosomesProbabilityIndex,
                        std::vector<std::vector<Key>> &intermediatePopulation);
  inline void generatePopulation(Graph<Key, double> &graph,
                                 int chromosomeNumber);
  inline void printBestSolution(Graph<Key, double> &graph, std::vector<std::pair<double, int>> &chromosomesScoreIndex);
  inline void printPopulation();
 public:
  explicit TSPSequential(Graph<Key, double> &graph);
  void Run(
         int chromosomeNumber,
           int generationNumber,
           double mutationRate,
           double crossoverRate,
           int workers,
           int seed) override;
};

template<typename Key>
TSPSequential<Key>::TSPSequential(Graph<Key, double> &graph): graph(graph) {}

template<typename Key>
void TSPSequential<Key>::Run(int chromosomeNumber,
                             int generationNumber,
                             double mutationRate,
                             double crossoverRate,
                             int workers,
                             int seed) {
  if (seed) {
    gen.seed(seed);
  }

  population.clear();
  population.resize(chromosomeNumber);
  auto start = std::chrono::system_clock::now();

  //! generate initial population
  generatePopulation(graph, chromosomeNumber);

  //! pairs of <fitness score, chromosomeIndex>
  std::vector<std::pair<double, int>> chromosomesScoreIndex(chromosomeNumber);
  std::vector<std::vector<Key>> intermediatePopulation(chromosomeNumber);
  double evaluationsAverage = 0;

  //! start genetic algorithm
  for (int generation = 1; generation <= generationNumber; generation++) {
    //std::cout << "Generation: " << generation << std::endl;
    evaluate(graph, chromosomesScoreIndex, evaluationsAverage);
    fitness(chromosomesScoreIndex, evaluationsAverage);
    printBestSolution(graph, chromosomesScoreIndex);
    selection(chromosomesScoreIndex, intermediatePopulation);
    crossover(intermediatePopulation, crossoverRate);
    mutation(intermediatePopulation, mutationRate);

    //! prepare for next generation
    evaluationsAverage = 0;
    intermediatePopulation.swap(population);
    intermediatePopulation.clear();
    chromosomesScoreIndex.clear();
    chromosomesScoreIndex.resize(chromosomeNumber);
    intermediatePopulation.resize(chromosomeNumber);
  }
  auto end = std::chrono::system_clock::now();
  std::cout << "Sequential time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms"
            << std::endl;
}

/*** Generate initial population. This is done by first creating vectors of pair<Key, probability>.
 *   Then all these vectors of pairs are sorted based on the random probability spawned to obtain
 *   the chromosome. In the end the population data structure is filled with copying the vector
 *   by taking only the Key element using std::transform
 *
 * @param initialPopulation Reference to the population data structure. Passed empty and later filled with chromosomes
 */
template<typename Key>
void TSPSequential<Key>::generatePopulation(Graph<Key, double> &graph,
                                            int chromosomeNumber) {
  //auto start = std::chrono::system_clock::now();
  std::vector<std::vector<std::pair<Key, double>>>
      populationToSort(chromosomeNumber);

  int nodesNumber = graph.GetNodesNumber();
  //! generate first population
  std::for_each(populationToSort.begin(),
                populationToSort.end(),
                [&](std::vector<std::pair<Key, double>> &chromosome) {
                  chromosome.resize(nodesNumber);
                  auto i = graph.GetMapIterator();

                  //! create chromosome probabilities
                  std::for_each(chromosome.begin(),
                                chromosome.end(),
                                [&](std::pair<Key, double> &cell) {
                                  cell.first = i->first;
                                  cell.second = unif(gen);
                                  i++;
                                });
                });
  /*
   //! print populationToSort for test purpose
   std::for_each(populationToSort.begin(),
                 populationToSort.end(),
                 [&](std::vector<std::pair<Key, double>> &chromosome) {

                   std::for_each(chromosome.begin(),
                                 chromosome.end(),
                                 [&](std::pair<Key, double> &cell) {
                                   std::cout << "Key: " << cell.first << " Value: " << cell.second << std::endl;
                                 });
                   std::cout << std::endl;
                 });
   */

  //! Iterate over the chromosomes and sorting in descending order (ordine decrescente) over to the probability
  std::for_each(populationToSort.begin(),
                populationToSort.end(),
                [&](std::vector<std::pair<Key, double>> &chromosome) {
                  std::sort(chromosome.begin(),
                            chromosome.end(),
                            [&](const std::pair<Key, double> &a, const std::pair<Key, double> &b) {
                              return a.second > b.second;
                            });
                });
  /*
  //! print populationToSort for test purpose
  std::for_each(populationToSort.begin(),
                populationToSort.end(),
                [&](std::vector<std::pair<Key, double>> &chromosome) {

                  std::for_each(chromosome.begin(),
                                chromosome.end(),
                                [&](std::pair<Key, double> &cell) {
                                  std::cout << "Key: " << cell.first << " Value: " << cell.second << std::endl;
                                });
                  std::cout << std::endl;
                }); */

  //! generate first population by taking Key values of each chromosomes generated
  auto populationToSortIt = populationToSort.begin();
  std::for_each(population.begin(),
                population.end(),
                [&](std::vector<Key> &chromosome) {
                  std::transform(populationToSortIt->begin(),
                                 populationToSortIt->end(),
                                 std::back_inserter(chromosome),
                                 [](const std::pair<Key, double> &pair) {
                                   return pair.first;
                                 });
                  populationToSortIt++;
                });
  //auto end = std::chrono::system_clock::now();
  //std::cout << "Generation time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

}

/*** For each chromosome compute the fitness score
 *   by summing all the edges value
 */
template<typename Key>
void TSPSequential<Key>::evaluate(Graph<Key, double> &graph,
                                  std::vector<std::pair<double, int>> &chromosomesScoreIndex,
                                  double &evaluationsAverage) {
  //auto start = std::chrono::system_clock::now();
  int chromosomeSize = population.begin()->size();
  auto chromosomesScoreIndexIt = chromosomesScoreIndex.begin();
  int iteratorIndex = 0;
  std::for_each(population.begin(),
                population.end(),
                [&](std::vector<Key> &chromosome) {
                  double chromosomeScore = 0;
                  for (int i = 0; i < chromosomeSize - 1; i++) {
                    chromosomeScore += graph.GetEdgeValue(chromosome[i], chromosome[i + 1]);
                  }
                  //! retrieve last edge value and store fitness score
                  chromosomesScoreIndexIt->first += chromosomeScore
                      + graph.GetEdgeValue(chromosome[0], chromosome[chromosomeSize - 1]);
                  chromosomesScoreIndexIt->second = iteratorIndex;
                  chromosomesScoreIndex[iteratorIndex].second = iteratorIndex;

                  evaluationsAverage += chromosomesScoreIndex[iteratorIndex].first;
                  chromosomesScoreIndexIt++;
                  iteratorIndex++;
                });
  evaluationsAverage /= chromosomesScoreIndex.size();
  //auto end = std::chrono::system_clock::now();
  //std::cout << "Evaluate time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

}

/*** Compute reproductive opportunity for each pair of chromosomesScoreIndex.
 *   This is done dividing the score (in our case the length of the path)
 *   by the average score obtained
 */
template<typename Key>
void TSPSequential<Key>::fitness(std::vector<std::pair<double, int>> &chromosomesScoreIndex,
                                 double evaluationsAverage) {
  /*
  //! print chromosome index with fitness value for test purpose
  std::for_each(chromosomesScoreIndex.begin(),
                chromosomesScoreIndex.end(), [&](std::pair<double, int> &currentIndex) {
        std::cout << "Evaluation: " << currentIndex.first << " Index: " << currentIndex.second << std::endl;
      }); */
  //auto start = std::chrono::system_clock::now();
  //printBestSolution(chromosomesScoreIndex); //! for test purpose
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

  //auto end = std::chrono::system_clock::now();
  //std::cout << "Fitness time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms"  << std::endl;

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
void TSPSequential<Key>::selection(std::vector<std::pair<double, int>> &chromosomesProbabilityIndex,
                                   std::vector<std::vector<Key>> &intermediatePopulation) {
  //auto start = std::chrono::system_clock::now();
  double probabilitiesSum = 0;
  std::vector<double> randomProbabilities(chromosomesProbabilityIndex.size());
  std::vector<int> indexNewPopulation(chromosomesProbabilityIndex.size());

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
  std::generate(randomProbabilities.begin(),
                randomProbabilities.end(),
                [&]() { return probabilitySpaceInterval(gen); });
  std::sort(randomProbabilities.begin(), randomProbabilities.end());

  //! fill the vector of new population indexes
  auto randomProbabilitiesIt = randomProbabilities.begin();
  auto indexNewPopulationIt = indexNewPopulation.begin();
  for (std::pair<double, int> &chromosome: chromosomesProbabilityIndex) {
    for (; randomProbabilitiesIt != randomProbabilities.end() && *randomProbabilitiesIt <= chromosome.first;
           randomProbabilitiesIt++, indexNewPopulationIt++) {
      *indexNewPopulationIt = chromosome.second;
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
  //auto end = std::chrono::system_clock::now();
  //std::cout << "Selection time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

}

/*** Implement crossover by picking up two random indexes in the list i and j with i<j and substituting the segment
 *   from i to j for each pair of chromosomes in the population
 */
template<typename Key>
void TSPSequential<Key>::crossover(std::vector<std::vector<Key>> &intermediatePopulation, const double crossoverRate) {
  //auto start = std::chrono::system_clock::now();
  std::random_shuffle(intermediatePopulation.begin(), intermediatePopulation.end());
  std::uniform_int_distribution<int> indexSpaceInterval(0, intermediatePopulation.begin()->size() - 1);
  std::vector<std::vector<Key>> crossoverPopulation(intermediatePopulation.size());

  //! swapping indexes sorted in increasing order
  std::pair<int, int> swappingIndexes(indexSpaceInterval(gen), indexSpaceInterval(gen));
  if (swappingIndexes.first > swappingIndexes.second) {
    std::swap(swappingIndexes.first, swappingIndexes.second);
  }

  //! apply crossover and save the results in crossoverPopulation
  auto intermediatePopulationNextIt = intermediatePopulation.begin() + 1;
  auto intermediatePopulationIt = intermediatePopulation.begin();
  auto crossoverPopulationIt = crossoverPopulation.begin();

  /*
  //! print intermediate population before crossover for test purpose
  std::for_each(intermediatePopulation.begin(),
                intermediatePopulation.end(),
                [&](std::vector<Key> &chromosome) {
                  std::for_each(chromosome.begin(),
                                chromosome.end(),
                                [&](Key &cell) {
                                  std::cout << cell << std::endl;
                                });
                  std::cout << std::endl;
                });
                */
  if (unif(gen) <= crossoverRate) {
    //! Apply Crossover between first and last element
    std::copy(intermediatePopulation.rbegin()->begin() + swappingIndexes.first,
              intermediatePopulation.rbegin()->begin() + swappingIndexes.second + 1,
              std::back_inserter(*crossoverPopulationIt));
    std::copy_if(intermediatePopulationIt->begin(),
                 intermediatePopulationIt->end(),
                 std::back_inserter(*crossoverPopulationIt),
                 [&](Key &currentKey) {
                   //! discard key if already present (iterator went through all crossoverPopulation array)
                   return
                       std::find(crossoverPopulationIt->begin(), crossoverPopulationIt->end(), currentKey)
                           == crossoverPopulationIt->end();
                 });
  } else {
    std::copy(intermediatePopulation.rbegin()->begin(),
              intermediatePopulation.rbegin()->end(),
              std::back_inserter(*crossoverPopulationIt));
  }
  crossoverPopulationIt++;

  //! Apply crossover to all the chromosomes (except last and first, already done)
  for (; intermediatePopulationNextIt != intermediatePopulation.end(); intermediatePopulationNextIt++) {
    if (unif(gen) <= crossoverRate) {
      std::copy(intermediatePopulationIt->begin() + swappingIndexes.first,
                intermediatePopulationIt->begin() + swappingIndexes.second + 1,
                std::back_inserter(*crossoverPopulationIt));
      std::copy_if(intermediatePopulationNextIt->begin(),
                   intermediatePopulationNextIt->end(),
                   std::back_inserter(*crossoverPopulationIt),
                   [&](Key &currentKey) {
                     //! discard key if already present (iterator went through all crossoverPopulation array)
                     return
                         std::find(crossoverPopulationIt->begin(), crossoverPopulationIt->end(), currentKey)
                             == crossoverPopulationIt->end();
                   });
    } else {
      std::copy(intermediatePopulationIt->begin(),
                intermediatePopulationIt->end(),
                std::back_inserter(*crossoverPopulationIt));
    }
    intermediatePopulationIt++;
    crossoverPopulationIt++;
  }
  crossoverPopulation.swap(intermediatePopulation);
  //auto end = std::chrono::system_clock::now();
  //std::cout << "Crossover time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

  /*
  //! print intermediate population after crossover for test purpose
  std::for_each(intermediatePopulation.begin(),
                intermediatePopulation.end(),
                [&](std::vector<Key> &chromosome) {
                  std::for_each(chromosome.begin(),
                                chromosome.end(),
                                [&](Key &cell) {
                                  std::cout << cell << std::endl;
                                });
                  std::cout << std::endl;
                });*/
  /*
 //! iterate over the intermediate population and do crossover
 std::for_each(swappingIndexes.begin(),
               swappingIndexes.end(), [&](std::pair<int, int> &swappingIndexes) {
           std::cout << "First: " << swappingIndexes.first << " Last: " << swappingIndexes.second << std::endl;
     }); */
}

/*** Mutate the chromosome by swapping their inner components according to mutation rate
 *
 * @param intermediatePopulation
 */
template<typename Key>
void TSPSequential<Key>::mutation(std::vector<std::vector<Key>> &intermediatePopulation, const double mutationRate) {
  //auto start = std::chrono::system_clock::now();
  std::uniform_int_distribution<int> randomIndexes(0, intermediatePopulation.begin()->size() - 1);

  std::for_each(intermediatePopulation.begin(),
                intermediatePopulation.end(),
                [&](std::vector<Key> &chromosome) {
                  /*
                  //! print before mutation for test purpose
                  std::for_each(chromosome.begin(),
                                chromosome.end(), [&](std::pair<Key, double> &currentIndex) {
                        std::cout << "Index: " << currentIndex.first << " Probability: " << currentIndex.second
                                  << std::endl;
                      }); */
                  //!  swap two random indexes of the chromosomes
                  if (unif(gen) <= mutationRate) {
                    std::swap(chromosome[randomIndexes(gen)], chromosome[randomIndexes(gen)]);
                  }
                  /*
                  //! print after mutation for test purpose
                  std::for_each(chromosome.begin(),
                                chromosome.end(), [&](std::pair<Key, double> &currentIndex) {
                        std::cout << "Index: " << currentIndex.first << " Probability: " << currentIndex.second
                                  << std::endl;
                      });
                      */
                });

  //auto end = std::chrono::system_clock::now();
  //std::cout << "Mutation time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

}

/*** Print population for test purpose
 *
 */
template<typename Key>
void TSPSequential<Key>::printPopulation() {
  //! print current population
  std::cout << "Current population " << std::endl;
  std::for_each(population.begin(),
                population.end(),
                [&](std::vector<Key> &chromosome) {
                  std::for_each(chromosome.begin(),
                                chromosome.end(),
                                [&](Key &cell) {
                                  std::cout << cell << std::endl;
                                });
                  std::cout << std::endl;
                });
}

/** Print current best solution
 *
 */
template<typename Key>
void TSPSequential<Key>::printBestSolution(Graph<Key, double> &graph,
                                           std::vector<std::pair<double, int>> &chromosomesScoreIndex) {
  std::vector<Key>& currentBestSolution = population.at(chromosomesScoreIndex.rbegin()->second);
  double chromosomeScore = 0;
  for (int i = 0; i < currentBestSolution.size() - 1; i++) {
    chromosomeScore += graph.GetEdgeValue(currentBestSolution[i], currentBestSolution[i + 1]);
  }
  chromosomeScore += graph.GetEdgeValue(currentBestSolution[0], currentBestSolution[currentBestSolution.size() - 1]);
  std::cout << chromosomeScore << std::endl;
}
#endif //GENETICTSP_SRC_TSP_SEQUENTIALTSP_H_