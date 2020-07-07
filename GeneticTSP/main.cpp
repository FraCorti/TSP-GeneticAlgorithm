#include "src/graph/graph.h"
#include "src/tsp/sequentialTSP.h"
#include "src/tsp/parallelTSP.h"
#include "src/tsp/fastflowTSP.h"

int main(int argc, char *argv[]) {
  if (argc < 7) {
    std::cout << "Usage is " << argv[0]
              << " nodesNumber chromosomesNumber generationNumber mutationRate crossoverRate workersNumber [seed]"
              << std::endl;
    return (-1);
  }

  int nodesNumber = std::atoi(argv[1]);
  int chromosomesNumber = std::atoi(argv[2]);
  int generationNumber = std::atoi(argv[3]);
  double mutationRate = std::atof(argv[4]);;
  double crossoverRate = std::atof(argv[5]);;
  int workerNumber = std::atoi(argv[6]);
  int seed = 0;
  std::string filePath;

  if (argv[7]) {
    seed = std::atoi(argv[7]);
  }
  if (argv[8]) {
    filePath = argv[8];
  }

  Graph<> graph(nodesNumber, seed);

  TSPSequential<int, double> tspSequential(graph);
  tspSequential.Run(chromosomesNumber, generationNumber, mutationRate, crossoverRate, workerNumber, seed);

  TSPParallel<int, double> tspParallel(graph);
  tspParallel.Run(chromosomesNumber, generationNumber, mutationRate, crossoverRate, workerNumber, seed);

  TSPFastflow<int, double> tspFastflow(graph);
  tspFastflow.Run(chromosomesNumber, generationNumber, mutationRate, crossoverRate, workerNumber, seed);
  return 0;
}
