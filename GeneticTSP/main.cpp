#include "src/graph/graph.h"
#include "src/tsp/sequentialTSP.h"

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cout << "Usage is " << argv[0]
              << " nodesNumber nw filePath" << std::endl;
    return (-1);
  }

  int totalNodesNumber = std::atoi(argv[1]);
  int workerNumber = std::atoi(argv[2]);
  std::string filePath;

  if (argv[3]) {
    filePath = argv[3];
  }
  Graph<> graph(totalNodesNumber);
  TSPSequential<>tsp(graph);
  tsp.Run(5000, 5000, 0.1, 0.1, 1, 0);
  return 0;
}