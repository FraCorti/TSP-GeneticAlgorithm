#include "src/graph/graph.h"
#include "src/tsp/tsp.h"

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cout << "Usage is " << argv[0]
              << "nnodes nw filePath" << std::endl;
    return (-1);
  }

  int totalNodesNumber = std::atoi(argv[1]);
  int workerNumber = std::atoi(argv[2]);
  std::string filePath;

  if(argv[3]){
    filePath = argv[3];
  }
  Graph graph(totalNodesNumber);
  TSP tsp(graph, 1, 1);
  return 0;
}