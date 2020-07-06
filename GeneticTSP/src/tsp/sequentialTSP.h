//
// Created by francesco on 31/05/20.
//

#ifndef GENETICTSP_SRC_TSP_SEQUENTIALTSP_H_
#define GENETICTSP_SRC_TSP_SEQUENTIALTSP_H_

#include "geneticAlgorithm.h"

template<typename Key = int, typename Value = double>
class TSPSequential : public GeneticAlgorithm {
 private:
 //! vector of chromosomes
  std::vector<std::vector<Key>> population;
  std::random_device rd;    // Will be used to obtain a seed for the random number engine
  std::mt19937 gen{rd()};   //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<double> unif{0, 1};
  Graph<Key, Value> &graph;

  void evaluate(Graph<Key, Value> &graph,
                std::vector<std::pair<precision , int>> &chromosomesScoreIndex,
                double &evaluationsAverage);
  void mutation(std::vector<std::vector<Key>> &intermediatePopulation, double mutationRate);
  void crossover(std::vector<std::vector<Key>> &intermediatePopulation, double crossoverRate);
  void fitness(std::vector<std::pair<precision, int>> &chromosomesScoreIndex, double evaluationsAverage) const;
  void selection(std::vector<std::pair<precision, int>> &chromosomesProbabilityIndex,
                 std::vector<std::vector<Key>> &intermediatePopulation);
  void generatePopulation(Graph<Key, Value> &graph,
                          int chromosomeNumber);
  void printBestSolution(Graph<Key, Value> &graph,
                         std::vector<std::pair<precision, int>> &chromosomesScoreIndex);
  void printPopulation() const;
 public:
  explicit TSPSequential(Graph<Key, Value> &graph);
  void Run(
      int chromosomeNumber,
      int generationNumber,
      double mutationRate,
      double crossoverRate,
      int workers,
      int seed) override;
};

template<typename Key, typename Value>
TSPSequential<Key, Value>::TSPSequential(Graph<Key, Value> &graph): graph(graph) {}

template<typename Key, typename Value>
void TSPSequential<Key, Value>::Run(const int chromosomeNumber,
                                    const int generationNumber,
                                    const double mutationRate,
                                    const double crossoverRate,
                                    const int workers,
                                    const int seed) {
  if (seed) {
    gen.seed(seed);
  }

  population.clear();
  population.resize(chromosomeNumber);
  auto start = std::chrono::system_clock::now();

  //! generate initial population
  generatePopulation(graph, chromosomeNumber);

  //! pairs of <fitness score, chromosomeIndex>
  std::vector<std::pair<precision, int>> chromosomesScoreIndex(chromosomeNumber);
  std::vector<std::vector<Key>> intermediatePopulation(chromosomeNumber);
  double evaluationsAverage = 0;

  //! start genetic algorithm
  for (size_t generation = 1; generation <= generationNumber; generation++) {
    evaluate(graph, chromosomesScoreIndex, evaluationsAverage);

    fitness(chromosomesScoreIndex, evaluationsAverage);
    //printBestSolution(graph, chromosomesScoreIndex); //! for test purpose
    selection(chromosomesScoreIndex, intermediatePopulation);

    crossover(intermediatePopulation, crossoverRate);
    mutation(intermediatePopulation, mutationRate);

    //! setup next generation
    evaluationsAverage = 0;
    population.swap(intermediatePopulation);
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
template<typename Key, typename Value>
void TSPSequential<Key, Value>::generatePopulation(Graph<Key, Value> &graph,
                                                   const int chromosomeNumber) {
  //auto start = std::chrono::system_clock::now();
  std::vector<std::vector<std::pair<Key, Value>>>
      populationToSort(chromosomeNumber);

  int nodesNumber = graph.GetNodesNumber();

  //! generate first population
  std::for_each(populationToSort.begin(),
                populationToSort.end(),
                [&](std::vector<std::pair<Key, Value>> &chromosome) -> void {
                  chromosome.resize(nodesNumber);
                  auto i = graph.GetMapIterator();

                  //! create chromosome probabilities
                  std::for_each(chromosome.begin(),
                                chromosome.end(),
                                [&](std::pair<Key, Value> &cell) -> void {
                                  cell.first = i->first;
                                  cell.second = unif(gen);
                                  i++;
                                });

                  //! sorting the chromosome in descending order over the probability (move out if you want to parallelize)
                  std::sort(chromosome.begin(),
                            chromosome.end(),
                            [](const std::pair<Key, Value> &a, const std::pair<Key, Value> &b) -> bool {
                              return a.second > b.second;
                            });

                });

  //! generate first population by taking Key values of each chromosomes generated
  auto populationToSortIt = populationToSort.begin();
  std::for_each(population.begin(),
                population.end(),
                [&populationToSortIt](std::vector<Key> &chromosome) {
                  std::transform(populationToSortIt->begin(),
                                 populationToSortIt->end(),
                                 std::back_inserter(chromosome),
                                 [](const std::pair<Key, Value> &pair) -> Key {
                                   return pair.first;
                                 });
                  populationToSortIt++;
                });
  auto end = std::chrono::system_clock::now();
  //std::cout << "Generation sequential time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
  //printPopulation();

}

/*** For each chromosome compute the fitness score
 *   by summing all the edges value
 */
template<typename Key, typename Value>
void TSPSequential<Key, Value>::evaluate(Graph<Key, Value> &graph,
                                         std::vector<std::pair<precision , int>> &chromosomesScoreIndex,
                                         double &evaluationsAverage) {
  //auto start = std::chrono::system_clock::now();
  int chromosomeSize = population.begin()->size();
  auto chromosomesScoreIndexIt = chromosomesScoreIndex.begin();
  int currentChromosomePosition = 0;
  std::for_each(population.begin(),
                population.end(),
                [&](const std::vector<Key> &chromosome) {
                  double chromosomeScore = 0;
                  for (int i = 0; i < chromosomeSize - 1; i++) {
                    chromosomeScore += graph.GetEdgeValue(chromosome.at(i), chromosome.at(i + 1));
                  }
                  //! retrieve edge value between last and first node and store fitness score
                  chromosomesScoreIndexIt->first = chromosomeScore
                      + graph.GetEdgeValue(chromosome.at(0), chromosome.at(chromosomeSize - 1));
                  chromosomesScoreIndexIt->second = currentChromosomePosition;

                  evaluationsAverage += chromosomesScoreIndexIt->first;
                  chromosomesScoreIndexIt++;
                  currentChromosomePosition++;
                });
  evaluationsAverage /= chromosomesScoreIndex.size();
  //auto end = std::chrono::system_clock::now();
  //std::cout << "Evaluate sequential time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

}

/*** Compute reproductive opportunity for each pair of chromosomesScoreIndex.
 *   This is done dividing the score (in our case the length of the path)
 *   by the average score obtained during the current generation
 */
template<typename Key, typename Value>
void TSPSequential<Key, Value>::fitness(std::vector<std::pair<precision, int>> &chromosomesScoreIndex,
                                        const double evaluationsAverage) const {
  //auto start = std::chrono::system_clock::now();
  //! for each chromosome compute and store probability of survive
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

  //auto end = std::chrono::system_clock::now();
  //std::cout << "Fitness time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms"  << std::endl;
}

/*** Produce the intermediate population based on the probability of reproduction computed before then
 *   the new intermediate population is saved inside population
 *
 *   @param chromosomesProbabilityIndex Vector of <reproduceProbability, chromosomeIndex> sorted in increasing order
 */
template<typename Key, typename Value>
void TSPSequential<Key, Value>::selection(std::vector<std::pair<precision, int>> &chromosomesProbabilityIndex,
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

  //! generate random numbers inside probability space and sort them
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
  //auto end = std::chrono::system_clock::now();
  //std::cout << "Selection time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

}

/*** Implement crossover by picking up two random indexes in the list i and j with i<j and substituting the segment
 *   from i to j for each pair of chromosomes in the population
 */
template<typename Key, typename Value>
void TSPSequential<Key, Value>::crossover(std::vector<std::vector<Key>> &intermediatePopulation,
                                          const double crossoverRate) {
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

  if (unif(gen) <= crossoverRate) {
    //! Apply Crossover between first and last element
    std::copy(intermediatePopulation.rbegin()->begin() + swappingIndexes.first,
              intermediatePopulation.rbegin()->begin() + swappingIndexes.second + 1,
              std::back_inserter(*crossoverPopulationIt));
    std::copy_if(intermediatePopulationIt->begin(),
                 intermediatePopulationIt->end(),
                 std::back_inserter(*crossoverPopulationIt),
                 [&](Key &currentKey) -> bool {
                   //! discard key if already present
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

  //! Apply crossover between current chromosome (intermediatePopulationIt) and next one (intermediatePopulationNextIt)
  for (; intermediatePopulationNextIt != intermediatePopulation.end(); intermediatePopulationNextIt++) {
    if (unif(gen) <= crossoverRate) {
      std::copy(intermediatePopulationIt->begin() + swappingIndexes.first,
                intermediatePopulationIt->begin() + swappingIndexes.second + 1,
                std::back_inserter(*crossoverPopulationIt));
      std::copy_if(intermediatePopulationNextIt->begin(),
                   intermediatePopulationNextIt->end(),
                   std::back_inserter(*crossoverPopulationIt),
                   [&](Key &currentKey) -> bool {
                     //! discard key if already present
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

}

/*** Mutate the chromosome by swapping their inner components according to mutation rate
 *
 * @param intermediatePopulation
 */
template<typename Key, typename Value>
void TSPSequential<Key, Value>::mutation(std::vector<std::vector<Key>> &intermediatePopulation,
                                         const double mutationRate) {
  //auto start = std::chrono::system_clock::now();
  std::uniform_int_distribution<int> randomIndexes(0, intermediatePopulation.begin()->size() - 1);

  std::for_each(intermediatePopulation.begin(),
                intermediatePopulation.end(),
                [&](std::vector<Key> &chromosome) -> void {
                  //!  swap two random indexes of the chromosomes
                  if (unif(gen) <= mutationRate) {
                    std::swap(chromosome.at(randomIndexes(gen)), chromosome.at(randomIndexes(gen)));
                  }
                });

  auto end = std::chrono::system_clock::now();
  //std::cout << "Mutation time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

}

/** Print current best solution
 *
 */
template<typename Key, typename Value>
void TSPSequential<Key, Value>::printBestSolution(Graph<Key, Value> &graph,
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

/*** Print population for test purpose
 *
 */
template<typename Key, typename Value>
void TSPSequential<Key, Value>::printPopulation() const {
  //! print current population
  std::cout << "Current sequential population " << std::endl;
  std::for_each(population.begin(),
                population.end(),
                [](const std::vector<Key> &chromosome) -> void {
                  std::for_each(chromosome.begin(),
                                chromosome.end(),
                                [&](const Key &cell) {
                                  std::cout << cell << std::endl;
                                });
                  std::cout << std::endl;
                });
}
#endif //GENETICTSP_SRC_TSP_SEQUENTIALTSP_H_