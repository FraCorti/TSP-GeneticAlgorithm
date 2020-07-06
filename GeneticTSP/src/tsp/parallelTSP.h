//
// Created by francesco on 27/06/20.
//

#ifndef GENETICTSP_SRC_TSP_PARALLELTSP_H_
#define GENETICTSP_SRC_TSP_PARALLELTSP_H_

#include <thread>
#include <future>
#include "geneticAlgorithm.h"

template<typename Key = int, typename Value = double>
class TSPParallel : public GeneticAlgorithm {
 private:

  //! vector of chromosomes
  std::vector<std::vector<Key>> population;
  std::vector<std::future<void>> workers;
  std::vector<std::pair<int, int>> chunks;

  std::random_device rd;    // Will be used to obtain a seed for the random number engine
  std::mt19937 gen{rd()};   //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<double> unif{0, 1};
  Graph<Key, Value> &graph;

  void parallelEvaluate(Graph<Key, Value> &graph,
                        std::vector<std::pair<precision, int>> &chromosomesScoreIndex,
                        double &evaluationsAverage);
  void mutation(std::vector<std::vector<Key>> &intermediatePopulation, double mutationRate);
  void parallelCrossover(std::vector<std::vector<Key>> &intermediatePopulation, double crossoverRate);
  void fitness(std::vector<std::pair<precision, int>> &chromosomesScoreIndex, double evaluationsAverage) const;
  void selection(std::vector<std::pair<precision, int>> &chromosomesProbabilityIndex,
                 std::vector<std::vector<Key>> &intermediatePopulation);
  void parallelGeneratePopulation(Graph<Key, Value> &graph,
                                  int chromosomeNumber);
  void printBestSolution(Graph<Key, Value> &graph,
                         std::vector<std::pair<precision, int>> &chromosomesScoreIndex);
  void setupComputation(int workersNumber, int chromosomeNumber);
  void printPopulation(const std::vector<std::vector<Key>> &population_) const;
 public:
  explicit TSPParallel(Graph<Key, Value> &graph);
  void Run(
      int chromosomeNumber,
      int generationNumber,
      double mutationRate,
      double crossoverRate,
      int workers_,
      int seed) override;

};

template<typename Key, typename Value>
void TSPParallel<Key, Value>::Run(const int chromosomeNumber,
                                  const int generationNumber,
                                  const double mutationRate,
                                  const double crossoverRate,
                                  const int workers_,
                                  const int seed) {
  if (seed) {
    gen.seed(seed);
  }

  population.clear();
  population.resize(chromosomeNumber);

  //! setup parallel computation
  setupComputation(workers_, chromosomeNumber);

  auto start = std::chrono::system_clock::now();

  //! generate initial population
  parallelGeneratePopulation(graph, chromosomeNumber);

  //! pairs of <fitness score, chromosomeIndex>
  std::vector<std::pair<precision, int>> chromosomesScoreIndex(chromosomeNumber);
  std::vector<std::vector<Key>> intermediatePopulation(chromosomeNumber);
  double evaluationsAverage = 0;

  //! start genetic algorithm
  for (size_t generation = 1; generation <= generationNumber; generation++) {
    parallelEvaluate(graph, chromosomesScoreIndex, evaluationsAverage);

    fitness(chromosomesScoreIndex, evaluationsAverage);
    //printBestSolution(graph, chromosomesScoreIndex);  // plot the convergence
    selection(chromosomesScoreIndex, intermediatePopulation);

    parallelCrossover(intermediatePopulation, crossoverRate);
    mutation(intermediatePopulation, mutationRate);

    //! setup next generation
    evaluationsAverage = 0;
    population.swap(intermediatePopulation);
  }
  auto end = std::chrono::system_clock::now();
  std::cout << "Standard C++ thread time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
  << "ms" << std::endl;
}

/*** Validate the number of worker passed by the user and allocate vector of futures
 *
 * @param workers Number of workers passed by the user
 */
template<typename Key, typename Value>
void TSPParallel<Key, Value>::setupComputation(const int workersNumber, const int chromosomeNumber) {
  int realWorkers =
      (workersNumber < std::thread::hardware_concurrency() - 1) ? workersNumber :
      std::thread::hardware_concurrency()
          - 1;

  if (realWorkers > chromosomeNumber) {
    realWorkers = chromosomeNumber;
  }

  workers.resize(realWorkers);
  chunks.resize(realWorkers);

  int workerCellsToCompute = std::floor(chromosomeNumber / realWorkers);
  int remainedElements = chromosomeNumber % realWorkers;
  int index = 0;

  //! setup chunk for workers
  std::for_each(chunks.begin(), chunks.end(), [&](std::pair<int, int> &chunk) {
    chunk.first = index;
    if (remainedElements) {
      chunk.second = chunk.first + workerCellsToCompute;
      remainedElements--;
    } else {
      chunk.second = chunk.first + workerCellsToCompute - 1;
    }
    index = chunk.second + 1;
  });
}

/*** Parallel generation of the population
 *
 */
template<typename Key, typename Value>
void TSPParallel<Key, Value>::parallelGeneratePopulation(Graph<Key, Value> &graph, int chromosomeNumber) {
  //auto start = std::chrono::system_clock::now();

  std::vector<std::vector<std::pair<Key, Value>>>
      populationToSort(chromosomeNumber);
  int nodesNumber = graph.GetNodesNumber();

  auto chunksIt = chunks.begin();

  std::for_each(populationToSort.begin(),
                populationToSort.end(),
                [&](std::vector<std::pair<Key, Value>> &chromosome) -> void {
                  chromosome.resize(nodesNumber);
                });

  std::for_each(workers.begin(), workers.end(), [&](std::future<void> &worker) {
    worker = std::async(std::launch::async, [&, start = chunksIt->first, end = chunksIt->second] {

      for (int chromosomeIndex = start; chromosomeIndex <= end; chromosomeIndex++) {
        std::vector<std::pair<Key, Value>> &chromosome = populationToSort.at(chromosomeIndex);
        auto mapIt = graph.GetMapIterator();

        //! create chromosome probabilities
        std::for_each(chromosome.begin(),
                      chromosome.end(),
                      [&](std::pair<Key, Value> &cell) -> void {
                        cell.first = mapIt->first;
                        cell.second = unif(gen);
                        mapIt++;
                      });

        //! sorting chromosome in descending order over the probability
        std::sort(chromosome.begin(),
                  chromosome.end(),
                  [](const std::pair<Key, Value> &a, const std::pair<Key, Value> &b) -> bool {
                    return a.second > b.second;
                  });

        //! fill population data structure with new chromosomes obtained after sorting
        std::transform(chromosome.begin(),
                       chromosome.end(),
                       std::back_inserter(population.at(chromosomeIndex)),
                       [](const std::pair<Key, Value> &pair) -> Key {
                         return pair.first;
                       });
      }
    });
    //! iterate over chunks data structure to assign indexes to the workers
    chunksIt++;
  });

  std::for_each(workers.begin(), workers.end(), [&](std::future<void> &worker) { worker.get(); });
  //auto end = std::chrono::system_clock::now();
  //std::cout << "Generation parallel time: "
  //        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
}

/*** Parallel implementation of evaluate operation

 */
template<typename Key, typename Value>
void TSPParallel<Key, Value>::parallelEvaluate(Graph<Key, Value> &graph,
                                               std::vector<std::pair<precision, int>> &chromosomesScoreIndex,
                                               double &evaluationsAverage) {
  //auto start = std::chrono::system_clock::now();
  auto chunksIt = chunks.begin();
  int chromosomeSize = population.begin()->size();
  std::for_each(workers.begin(), workers.end(), [&](std::future<void> &worker) {
    worker = std::async(std::launch::async,
                        [&, start = chunksIt->first, end = chunksIt->second, chromosomeSize = chromosomeSize] {
                          for (int chromosomeIndex = start; chromosomeIndex <= end; chromosomeIndex++) {
                            double chromosomeScore = 0;
                            std::vector<Key> &chromosome = population.at(chromosomeIndex);
                            for (int i = 0; i < chromosomeSize - 1; i++) {
                              chromosomeScore += graph.GetEdgeValue(chromosome.at(i), chromosome.at(i + 1));
                            }
                            chromosomesScoreIndex.at(chromosomeIndex).first = chromosomeScore
                                + graph.GetEdgeValue(chromosome.at(0), chromosome.at(chromosomeSize - 1));
                            chromosomesScoreIndex.at(chromosomeIndex).second = chromosomeIndex;
                          }
                        });
    //! iterate over chunks data structure to assign indexes to the workers
    chunksIt++;
  });
  std::for_each(workers.begin(), workers.end(), [&](std::future<void> &worker) { worker.get(); });
  double totalScore = 0;
  std::for_each(chromosomesScoreIndex.begin(),
                chromosomesScoreIndex.end(),
                [&totalScore](const std::pair<precision, int> scoreIndex) {
                  totalScore += scoreIndex.first;
                });
  evaluationsAverage = totalScore / chromosomesScoreIndex.size();
}

/***
 *
 * @param chromosomesScoreIndex
 * @param evaluationsAverage
 */
template<typename Key, typename Value>
void TSPParallel<Key, Value>::fitness(std::vector<std::pair<precision, int>> &chromosomesScoreIndex,
                                      double evaluationsAverage) const {
  std::for_each(chromosomesScoreIndex.begin(),
                chromosomesScoreIndex.end(),
                [&evaluationsAverage](std::pair<double, int> &chromosome) -> void {
                  //! formula adapted from "A genetic algorithm tutorial" Darrel Whitley
                  chromosome.first = evaluationsAverage / chromosome.first;
                });

  //! sort obtained probability in increasing order (ordine crescente)
  std::sort(chromosomesScoreIndex.begin(),
            chromosomesScoreIndex.end(),
            [](const std::pair<precision, int> &a, const std::pair<precision, int> &b) -> bool {
              return a.first < b.first;
            });
}

/***
 *
 * @param chromosomesProbabilityIndex
 * @param intermediatePopulation
 */
template<typename Key, typename Value>
void TSPParallel<Key, Value>::selection(std::vector<std::pair<precision, int>> &chromosomesProbabilityIndex,
                                        std::vector<std::vector<Key>> &intermediatePopulation) {
  //auto start = std::chrono::system_clock::now();
  double probabilitiesSum = 0;
  std::vector<precision> randomProbabilities(chromosomesProbabilityIndex.size());
  std::vector<int> indexNewPopulation(chromosomesProbabilityIndex.size());

  //! divide the probability space between the chromosome over their reproduction factor
  std::for_each(chromosomesProbabilityIndex.begin(),
                chromosomesProbabilityIndex.end(),
                [&probabilitiesSum](std::pair<precision, int> &chromosome) -> void {
                  probabilitiesSum += chromosome.first;
                  chromosome.first = probabilitiesSum;
                });

  //! generate random numbers inside probability space and sort them (to enforce data locality)
  std::uniform_real_distribution<double> probabilitySpaceInterval(0, probabilitiesSum);
  std::generate(randomProbabilities.begin(),
                randomProbabilities.end(),
                [&]() -> double { return probabilitySpaceInterval(gen); });
  std::sort(randomProbabilities.begin(), randomProbabilities.end());

  //! fill the vector of new population indexes
  auto randomProbabilitiesIt = randomProbabilities.begin();
  auto indexNewPopulationIt = indexNewPopulation.begin();
  for (std::pair<precision, int> &chromosome: chromosomesProbabilityIndex) {
    for (; randomProbabilitiesIt != randomProbabilities.end() && *randomProbabilitiesIt <= chromosome.first;
           randomProbabilitiesIt++, indexNewPopulationIt++) {
      *indexNewPopulationIt = chromosome.second;
    }
  }

  //! create intermediate population
  auto intermediatePopulationIt = intermediatePopulation.begin();
  std::for_each(indexNewPopulation.begin(),
                indexNewPopulation.end(), [&](int &index) -> void {
        *intermediatePopulationIt = population.at(index);
        intermediatePopulationIt++;
      });
}

/*** Parallel implementation (using std::async) of Crossover operation.
 *   Basically the same logic applied in the sequential version is used but here
 *   a chunk that defines the range of chromosome inside population data structure is passed.
 *
 *   @param intermediatePopulation
 *   @param crossoverRate
 */
template<typename Key, typename Value>
void TSPParallel<Key, Value>::parallelCrossover(std::vector<std::vector<Key>> &intermediatePopulation, double crossoverRate) {
  // auto start = std::chrono::system_clock::now();
  std::random_shuffle(intermediatePopulation.begin(), intermediatePopulation.end());
  std::uniform_int_distribution<int> indexSpaceInterval(0, intermediatePopulation.begin()->size() - 1);
  std::vector<std::vector<Key>> crossoverPopulation(intermediatePopulation.size());

  //! swapping indexes sorted in increasing order
  std::pair<int, int> swappingIndexes(indexSpaceInterval(gen), indexSpaceInterval(gen));
  if (swappingIndexes.first > swappingIndexes.second) {
    std::swap(swappingIndexes.first, swappingIndexes.second);
  }

  auto chunksIt = chunks.begin();
  std::for_each(workers.begin(), workers.end(), [&](std::future<void> &worker) {
    worker = std::async(std::launch::async,
                        [&, start = chunksIt->first, end = chunksIt->second, crossoverStart = swappingIndexes.first,
                            crossoverEnd = swappingIndexes.second] {

                          //! Apply Crossover between first and last element
                          if (unif(gen) <= crossoverRate) {
                            std::vector<Key> &firstChromosome = intermediatePopulation.at(start);
                            std::vector<Key> &lastChromosome = intermediatePopulation.at(end);
                            std::copy(lastChromosome.begin() + crossoverStart,
                                      lastChromosome.begin() + crossoverEnd + 1,
                                      std::back_inserter(crossoverPopulation.at(end)));
                            std::copy_if(firstChromosome.begin(),
                                         firstChromosome.end(),
                                         std::back_inserter(crossoverPopulation.at(end)),
                                         [&](Key &currentKey) -> bool {
                                           //! discard key if already present
                                           return
                                               std::find(crossoverPopulation.at(end).begin(),
                                                         crossoverPopulation.at(end).end(),
                                                         currentKey)
                                                   == crossoverPopulation.at(end).end();
                                         });
                          } else {
                            std::copy(intermediatePopulation.at(end).begin(),
                                      intermediatePopulation.at(end).end(),
                                      std::back_inserter(crossoverPopulation.at(end)));
                          }

                          //! Apply Crossover to inner chromosomes of the chunk
                          for (int chromosomeIndex = start; chromosomeIndex <= end - 1; chromosomeIndex++) {
                            std::vector<Key> &chromosome = intermediatePopulation.at(chromosomeIndex);
                            std::vector<Key> &nextChromosome = intermediatePopulation.at(chromosomeIndex + 1);
                            if (unif(gen) <= crossoverRate) {
                              std::copy(chromosome.begin() + crossoverStart,
                                        chromosome.begin() + crossoverEnd + 1,
                                        std::back_inserter(crossoverPopulation.at(chromosomeIndex)));
                              std::copy_if(nextChromosome.begin(),
                                           nextChromosome.end(),
                                           std::back_inserter(crossoverPopulation.at(chromosomeIndex)),
                                           [&](Key &currentKey) -> bool {
                                             //! discard key if already present
                                             return
                                                 std::find(crossoverPopulation.at(chromosomeIndex).begin(),
                                                           crossoverPopulation.at(chromosomeIndex).end(),
                                                           currentKey)
                                                     == crossoverPopulation.at(chromosomeIndex).end();
                                           });
                            } else {
                              std::copy(chromosome.begin(),
                                        chromosome.end(),
                                        std::back_inserter(crossoverPopulation.at(chromosomeIndex)));
                            }
                          }
                        });
    chunksIt++;
  });
  std::for_each(workers.begin(), workers.end(), [&](std::future<void> &worker) { worker.get(); });
  crossoverPopulation.swap(intermediatePopulation);
  //auto end = std::chrono::system_clock::now();
  //std::cout << "Crossover parallel time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
}

template<typename Key, typename Value>
void TSPParallel<Key, Value>::mutation(std::vector<std::vector<Key>> &intermediatePopulation, double mutationRate) {
  //auto start = std::chrono::system_clock::now();
  std::uniform_int_distribution<int> randomIndexes(0, intermediatePopulation.begin()->size() - 1);

  std::for_each(intermediatePopulation.begin(),
                intermediatePopulation.end(),
                [&](std::vector<Key> &chromosome) -> void {
                  //! swap two random cells
                  if (unif(gen) <= mutationRate) {
                    std::swap(chromosome.at(randomIndexes(gen)), chromosome.at(randomIndexes(gen)));
                  }
                });
  //auto end = std::chrono::system_clock::now();
  //std::cout << "Mutation parallel time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms"
  //<< std::endl;
}

/*** Useful for plot the convergence of the algorithm
 *
 */
template<typename Key, typename Value>
void TSPParallel<Key, Value>::printBestSolution(Graph<Key, Value> &graph,
                                                std::vector<std::pair<precision, int>> &chromosomesScoreIndex) {
  std::vector<Key> &currentBestSolution = population.at(chromosomesScoreIndex.rbegin()->second);
  double chromosomeScore = 0;
  for (int i = 0; i < currentBestSolution.size() - 1; i++) {
    chromosomeScore += graph.GetEdgeValue(currentBestSolution.at(i), currentBestSolution.at(i + 1));
  }
  chromosomeScore +=
      graph.GetEdgeValue(currentBestSolution.at(0), currentBestSolution.at(currentBestSolution.size() - 1));
  std::cout << chromosomeScore << std::endl;
}

/*** Print current population for test purpose
 *
 */
template<typename Key, typename Value>
void TSPParallel<Key, Value>::printPopulation(const std::vector<std::vector<Key>> &population_) const {
  std::for_each(population_.begin(),
                population_.end(),
                [](const std::vector<Key> &chromosome) -> void {
                  std::for_each(chromosome.begin(),
                                chromosome.end(),
                                [&](const Key &cell) {
                                  std::cout << cell << std::endl;
                                });
                  std::cout << std::endl;
                });
}
template<typename Key, typename Value>
TSPParallel<Key, Value>::TSPParallel(Graph<Key, Value>
                                     &graph): graph(graph) {}
#endif //GENETICTSP_SRC_TSP_PARALLELTSP_H_