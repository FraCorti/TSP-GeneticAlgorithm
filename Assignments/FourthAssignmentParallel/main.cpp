#include "ImageProcessing.h"


int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cout << "Usage is " << argv[0]
              << "directoryPath nw multiplyWork" << std::endl;
    return (-1);
  }

  //! directory containing the files
  std::string path = argv[1];;
  int workersNumber = atoi(argv[2]);
  int multiplyWork = atoi(argv[3]);

  ImageProcessing imageProcessing(path, workersNumber, multiplyWork);

  //! Sequential two stages work, used to benchmark the filters individually
  imageProcessing.SequentialTwoStage();

  //! Sequential one stage
  imageProcessing.SequentialOneStage();

  //! Parallel version OpenMP
  imageProcessing.ParallelOpenMP();

  //! Parallel version standard C++ threads
  imageProcessing.ParallelStandard();

  return 0;
}
