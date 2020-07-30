//
// Created by francesco on 13/06/20.
//

#ifndef GENETICTSP_SRC_TSP_FASTFLOWTSP_H_
#define GENETICTSP_SRC_TSP_FASTFLOWTSP_H_

#include "geneticAlgorithm.h"
#include "ff/ff.hpp"
#include <ff/pipeline.hpp>
#include <ff/farm.hpp>
#include <ff/pipeline.hpp>
#include <memory>
#include <thread>

template<typename Key = int, typename Value = double>
class TSPFastflow : public GeneticAlgorithm {
 private:
  struct ChromosomeGenerationWorker : ff::ff_node_t<std::pair<int, int>> {
   private:
    std::vector<std::vector<Key>> &population;
    std::vector<std::vector<std::pair<Key, Value>>> &initialGeneratedPopulation;
    Graph<Key, Value> &graph;
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<double> unif{0, 1};
   public:
    ChromosomeGenerationWorker(std::vector<std::vector<Key>> *population,
                               Graph<Key, Value> *graph,
                               std::vector<std::vector<std::pair<Key, Value>>> *initialGeneratedPopulation_) :
        population(*population), graph(*graph), initialGeneratedPopulation(*initialGeneratedPopulation_) {}
    std::pair<int, int> *svc(std::pair<int, int> *chunk) override {
      int start = chunk->first;
      int end = chunk->second;
      delete chunk;

      for (int chromosomeIndex = start; chromosomeIndex <= end; chromosomeIndex++) {
        std::vector<std::pair<Key, Value>> &chromosome = initialGeneratedPopulation.at(chromosomeIndex);
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
      this->ff_send_out(new std::pair<int, int>(start, end));
      return this->GO_ON;
    }
  };

  /*** Worker for crossover and mutation phase
   *
   */
  struct CrossoverMutationWorker : ff::ff_node_t<std::pair<int, int>> {
   private:
    const double crossoverRate;
    const double mutationRate;
    std::vector<std::vector<Key>> &population;
    std::vector<std::vector<Key>> &crossoverPopulation;
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<double> unif{0, 1};
    std::uniform_int_distribution<int> indexSpaceInterval;
   public:
    explicit CrossoverMutationWorker(const double crossover_rate,
                                     const double mutation_rate,
                                     std::vector<std::vector<Key>> *population,
                                     std::vector<std::vector<Key>> &crossover_population,
                                     const int graphNode)
        : crossoverRate(crossover_rate),
          mutationRate(mutation_rate),
          population(*population),
          crossoverPopulation(crossover_population), indexSpaceInterval(0, graphNode - 1) {
    }
    std::pair<int, int> *svc(std::pair<int, int> *chunk) override {
      int startIndex = chunk->first;
      int endIndex = chunk->second;

      //! generate swapping indexes
      int crossoverStart = indexSpaceInterval(gen);
      int crossoverEnd = indexSpaceInterval(gen);
      if (crossoverStart > crossoverEnd) {
        std::swap(crossoverStart, crossoverEnd);
      }

      //! Apply Crossover between first and last element of the chunk
      if (unif(gen) <= crossoverRate) {
        std::vector<Key> &firstChromosome = population.at(startIndex);
        std::vector<Key> &lastChromosome = population.at(endIndex);
        std::copy(lastChromosome.begin() + crossoverStart,
                  lastChromosome.begin() + crossoverEnd + 1,
                  std::back_inserter(crossoverPopulation.at(endIndex)));
        std::copy_if(firstChromosome.begin(),
                     firstChromosome.end(),
                     std::back_inserter(crossoverPopulation.at(endIndex)),
                     [&](Key &currentKey) -> bool {
                       //! discard key if already present
                       return
                           std::find(crossoverPopulation.at(endIndex).begin(),
                                     crossoverPopulation.at(endIndex).end(),
                                     currentKey)
                               == crossoverPopulation.at(endIndex).end();
                     });
      } else {
        std::copy(population.at(endIndex).begin(),
                  population.at(endIndex).end(),
                  std::back_inserter(crossoverPopulation.at(endIndex)));
      }

      //! Apply Mutation to the last chromosome of the chunk
      if (unif(gen) <= mutationRate) {
        std::swap(crossoverPopulation.at(endIndex).at(indexSpaceInterval(gen)),
                  crossoverPopulation.at(endIndex).at(indexSpaceInterval(gen)));
      }

      //! Apply Crossover and Mutation to the chromosomes of the chunk (except last one)
      for (int chromosomeIndex = startIndex; chromosomeIndex <= endIndex - 1; chromosomeIndex++) {
        std::vector<Key> &chromosome = population.at(chromosomeIndex);
        std::vector<Key> &nextChromosome = population.at(chromosomeIndex + 1);
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

        //! apply mutation to all the chromosomes of the chunk
        if (unif(gen) <= mutationRate) {
          std::swap(crossoverPopulation.at(chromosomeIndex).at(indexSpaceInterval(gen)),
                    crossoverPopulation.at(chromosomeIndex).at(indexSpaceInterval(gen)));
        }
      }

      //! store the results in population data structure
      for (int chromosomeIndex = startIndex; chromosomeIndex <= endIndex; chromosomeIndex++) {
        population.at(chromosomeIndex) = crossoverPopulation.at(chromosomeIndex);
      }

      this->ff_send_out(chunk);
      return this->GO_ON;
    }
  };

  //! emitter node for evaluate stage
  struct ChunksEmitter : ff::ff_node_t<std::pair<int, int>> {
   private:
    std::vector<std::pair<int, int>> &chunks;
    int generationsNumber;
   public:
    explicit ChunksEmitter(std::vector<std::pair<int, int>> *chunks, int generationsNumber_)
        : chunks(*chunks), generationsNumber(generationsNumber_) {}

    std::pair<int, int> *svc(std::pair<int, int> *currentGeneration) override {
      if (generationsNumber <= 0) {
        return this->EOS;
      }
      generationsNumber--;
      std::for_each(chunks.begin(), chunks.end(), [&](std::pair<int, int> chunk) {
        this->ff_send_out(new std::pair<int, int>(chunk.first, chunk.second));
      });
      return this->GO_ON;
    }
  };

   //! worker node for Evaluate stage
  struct EvaluateWorker : ff::ff_node_t<std::pair<int, int>> {
   private:
    std::vector<std::vector<Key>> &population;
    std::vector<std::pair<precision, int>> &chromosomesScoreIndex;
    Graph<Key, Value> &graph;
    const int chromosomeSize;
   public:
    explicit EvaluateWorker(std::vector<std::vector<Key>> *population,
                            std::vector<std::pair<precision, int>> &chromosomes_score_index,
                            Graph<Key, Value> *graph)
        : population(*population),
          chromosomesScoreIndex(chromosomes_score_index),
          graph(*graph),
          chromosomeSize(population->begin()->size()) {
    }

    std::pair<int, int> *svc(std::pair<int, int> *chunk) override {
      int start = chunk->first;
      int end = chunk->second;
      delete chunk;
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
      this->ff_send_out(new std::pair<int, int>(start, end));
      return this->GO_ON;
    }
  };

  struct CrossoverMutationCollector : ff::ff_node_t<std::pair<int, int>> {
   private:
    const int chunksNumber;
    const int chromosomeNumber;
    int currentChunksNumber;
    std::vector<std::vector<Key>> &crossoverPopulation;
   public:
    explicit CrossoverMutationCollector(std::vector<std::vector<Key>> &crossover_population,
                                        int chromosomeNumber,
                                        int chunks_number)
        : crossoverPopulation(crossover_population),
          chromosomeNumber(chromosomeNumber),
          chunksNumber(chunks_number),
          currentChunksNumber(chunks_number) {}
    std::pair<int, int> *svc(std::pair<int, int> *chunk) override {
      currentChunksNumber--;

      //! if all the chunks have been received then proceed to the next stage
      if (currentChunksNumber == 0) {
        crossoverPopulation.clear();
        crossoverPopulation.resize(chromosomeNumber);
        currentChunksNumber = chunksNumber;
        this->ff_send_out(chunk);
      }
      return this->GO_ON;
    }
  };

  struct ChunksCollector : ff::ff_node_t<std::pair<int, int>> {
   private:
    const int chunksNumber;
    int currentChunksNumber;
   public:
    explicit ChunksCollector(int chunks_number) : chunksNumber(chunks_number), currentChunksNumber(chunks_number) {}
    std::pair<int, int> *svc(std::pair<int, int> *chunk) override {
      currentChunksNumber--;

      //! if all the chunks have been computed then proceed to the next stage
      if (currentChunksNumber == 0) {
        currentChunksNumber = chunksNumber;
        this->ff_send_out(chunk);
      }
      return this->GO_ON;
    }
  };

/**  PRE: chromosomesScoreIndex has been filled with pairs of <fitnessScore, chromosomeIndex>
 * */
  struct FitnessSelection : ff::ff_node_t<std::pair<int, int>> {
   private:
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<double> unif{0, 1};
    std::vector<std::vector<Key>> &population;
    std::vector<std::pair<precision, int>> &chromosomesScoreIndex;
    std::vector<std::vector<Key>> intermediatePopulation;
   public:
    explicit FitnessSelection(std::vector<std::vector<Key>> *population_,
                              std::vector<std::pair<precision, int>> &chromosomes_score_index,
                              const int chunks_number = 0,
                              int seed = 0) : population(*population_),
                                              chromosomesScoreIndex(chromosomes_score_index) {
      intermediatePopulation.resize(population_->size());
      if (seed) {
        gen.seed(seed);
      }
    }
    std::pair<int, int> *svc(std::pair<int, int> *task) override {
      delete task;
      double totalScore = 0;
      std::for_each(chromosomesScoreIndex.begin(),
                    chromosomesScoreIndex.end(),
                    [&totalScore](const std::pair<precision, int> scoreIndex) {
                      totalScore += scoreIndex.first;
                    });
      double evaluationsAverage = totalScore / chromosomesScoreIndex.size();
      //printBestSolution(chromosomesScoreIndex); // for convergence test purpose
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

      double probabilitiesSum = 0;
      std::vector<precision> randomProbabilities(chromosomesScoreIndex.size());
      std::vector<int> indexNewPopulation(chromosomesScoreIndex.size());

      //! divide the probability space between the chromosome over their reproduction factor
      std::for_each(chromosomesScoreIndex.begin(),
                    chromosomesScoreIndex.end(),
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
      for (std::pair<precision, int> &chromosome: chromosomesScoreIndex) {
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
      population.swap(intermediatePopulation);

      //! proceed to next stage
      this->ff_send_out(new std::pair<int, int>(1, 1));
      return this->GO_ON;
    }
  };

  //! vector of chromosomes
  std::vector<std::vector<Key>> population;
  std::random_device rd;
  std::mt19937 gen{rd()};
  std::uniform_real_distribution<double> unif{0, 1};
  Graph<Key, Value> &graph;
  std::vector<std::pair<int, int>> chunks;

  void generatePopulationFastflow(Graph<Key, Value> &graph,
                                  int chromosomeNumber, int workers, int nodesNumber);
  void setupComputation(int &workersNumber, int chromosomeNumber);
  static void printBestSolution(std::vector<std::pair<precision, int>> &chromosomesScoreIndex);
  void printPopulation() const;
 public:
  explicit TSPFastflow(Graph<Key, Value> &graph);
  void Run(
      int chromosomeNumber,
      int generationNumber,
      double mutationRate,
      double crossoverRate,
      int workers,
      int seed) override;
};

template<typename Key, typename Value>
TSPFastflow<Key, Value>::TSPFastflow(Graph<Key, Value> &graph): graph(graph) {}

template<typename Key, typename Value>
void TSPFastflow<Key, Value>::Run(int chromosomeNumber,
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
  setupComputation(workers, chromosomeNumber);

  auto start = std::chrono::system_clock::now();

  //! generate initial population using Fastflow farm
  generatePopulationFastflow(graph, chromosomeNumber, workers, graph.GetNodesNumber());

  //! pairs of <fitness score, chromosomeIndex>
  std::vector<std::pair<precision, int>> chromosomesScoreIndex(chromosomeNumber);
  std::vector<std::vector<Key>> crossoverPopulation(chromosomeNumber);

  //! evaluate stage farm
  ChunksEmitter evaluateEmitter(&chunks, generationNumber);
  ChunksCollector evaluateCollector(chunks.size());
  std::vector<std::unique_ptr<ff::ff_node>> evaluateWorkers;
  for (int i = 0; i < workers; ++i) {
    evaluateWorkers.push_back(std::unique_ptr<EvaluateWorker>(new EvaluateWorker(&population,
                                                                                 chromosomesScoreIndex,
                                                                                 &graph)));
  }
  ff::ff_Farm<std::pair<int, int>> evaluateFarm(std::move(evaluateWorkers), evaluateEmitter, evaluateCollector);

  //! intermediate stage
  FitnessSelection fitnessSelectionStage(&population, chromosomesScoreIndex);

  //! crossover and mutation farm
  ChunksEmitter crossoverMutationEmitter(&chunks, generationNumber);
  CrossoverMutationCollector crossoverMutationCollector(crossoverPopulation, chromosomeNumber, chunks.size());
  std::vector<std::unique_ptr<ff::ff_node>> crossoverMutationWorkers;
  for (int i = 0; i < workers; ++i) {
    crossoverMutationWorkers.push_back(std::unique_ptr<CrossoverMutationWorker>(new CrossoverMutationWorker(
        crossoverRate,
        mutationRate,
        &population,
        crossoverPopulation,
        graph.GetNodesNumber())));;
  }
  ff::ff_Farm<std::pair<int, int>>
      crossoverMutationFarm(std::move(crossoverMutationWorkers), crossoverMutationEmitter, crossoverMutationCollector);

  ff::ff_Pipe<std::pair<int, int>> pipe(evaluateFarm, fitnessSelectionStage, crossoverMutationFarm);
  pipe.wrap_around();
  if (pipe.run_and_wait_end() < 0) {
    ff::error("running pipe");
    return;
  }
  auto end = std::chrono::system_clock::now();
  std::cout << "Fastflow time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms "
            << std::endl;
}

/*** Setup chunk data structure and checks number of workers passed in input
 *
 * @param workersNumber
 * @param chromosomeNumber
 */
template<typename Key, typename Value>
void TSPFastflow<Key, Value>::setupComputation(int &workersNumber, int chromosomeNumber) {
  workersNumber =
      (workersNumber < std::thread::hardware_concurrency() - 1) ? workersNumber :
      std::thread::hardware_concurrency()
          - 1;

  if (workersNumber > chromosomeNumber) {
    workersNumber = chromosomeNumber;
  }

  chunks.resize(workersNumber);

  int workerCellsToCompute = std::floor(chromosomeNumber / workersNumber);
  int remainedElements = chromosomeNumber % workersNumber;
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

/*** Generate initial population
 *
 * @param graph Current graph considered
 * @param chromosomeNumber
 */
template<typename Key, typename Value>
void TSPFastflow<Key, Value>::generatePopulationFastflow(Graph<Key, Value> &graph,
                                                         int chromosomeNumber,
                                                         int workers,
                                                         int nodesNumber) {

  std::vector<std::vector<std::pair<Key, Value>>> initialGeneratedPopulation(chromosomeNumber);
  std::for_each(initialGeneratedPopulation.begin(),
                initialGeneratedPopulation.end(),
                [&](std::vector<std::pair<Key, Value>> &chromosome) -> void {
                  chromosome.resize(nodesNumber);
                });

  //! chromosomes creation farm
  ChunksEmitter generateChromosomesEmitter(&chunks, 1);
  ChunksCollector generateChromosomesCollector(chunks.size());
  std::vector<std::unique_ptr<ff::ff_node>> generationWorkers;
  for (int i = 0; i < workers; i++) {
    generationWorkers.push_back(std::unique_ptr<ChromosomeGenerationWorker>(new ChromosomeGenerationWorker(&population,
                                                                                                           &graph,
                                                                                                           &initialGeneratedPopulation)));
  }
  ff::ff_Farm<std::pair<int, int>>
      creationFarm(std::move(generationWorkers), generateChromosomesEmitter, generateChromosomesCollector);
  creationFarm.wrap_around();
  if (creationFarm.run_and_wait_end() < 0) {
    return;
  }

}

/*** This method is useful if you need to plot the convergence of the algorithm
 *
 * @param graph
 * @param chromosomesScoreIndex
 */
template<typename Key, typename Value>
void TSPFastflow<Key, Value>::printBestSolution(std::vector<std::pair<precision, int>> &chromosomesScoreIndex) {

  //! sort in ascending order
  std::sort(chromosomesScoreIndex.begin(),
            chromosomesScoreIndex.end(),
            [](const std::pair<precision, int> a, const std::pair<precision, int> b) {
              return a.first < b.first;
            });

  std::cout << chromosomesScoreIndex.begin()->first << std::endl;
}

/** Print current population (useful for test purpose)
 */
template<typename Key, typename Value>
void TSPFastflow<Key, Value>::printPopulation() const {
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

#endif //GENETICTSP_SRC_TSP_FASTFLOWTSP_H_