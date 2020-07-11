# TSP with Genetic Algorithms
The purpose of the project was to solve [Travelling Salesman Problem](https://en.wikipedia.org/wiki/Travelling_salesman_problem), known as *TSP*, by using a [Genetic Algorithms](https://ibug.doc.ic.ac.uk/media/uploads/documents/courses/GeneticAlgorithm-tutorial.pdf) approach. In addition to the sequential version used as baseline, two parallel versions were developed, one using standard C++ threads and one using [Fastflow](https://github.com/fastflow/fastflow#fastflow-high-performance-parallel-patterns-and-building-blocks-in-c). More information about the project and the design of the parallel versions can be found in the [report](https://github.com/FraCorti/TSP-GeneticAlgorithm/blob/master/GeneticTSP/report.pdf).

To test the correcteness of the three version developed, I saved for a run of each version the current best route for each generation of the algorithm. In the following plots the length of the best path decrease as the algorithm go on with the generation, so convergence algorithms were obtained. 


<img src="https://github.com/FraCorti/TSP-GeneticAlgorithm/GeneticTSP/report/benchmark/convergence/fastflowConvergenceREADME.png" width="300"/> <img src="https://github.com/FraCorti/TSP-GeneticAlgorithm/GeneticTSP/report/benchmark/convergence/sequentialConvergenceREADME.png" width="300"/> <img src="https://github.com/FraCorti/TSP-GeneticAlgorithm/GeneticTSP/report/benchmark/convergence/standardConvergenceREADME.png" width="300"/>


## Running the project
The project can be compiled and executed using [Cmake 3.11](https://cmake.org/) or [GCC](https://gcc.gnu.org/) (remember to turn on -O3 optimization). 

If CMake is installed given:
`mkdir build && cd build ` 

In the end to build the project type: 
`cmake --build .` 

Otherwise, to use GCC type: 

`g++ -std=c++17 -O3 main.cpp -o GeneticTSP -pthread `

To run the program correctly, you need to pass to the executable the following parameters: 

```./GeneticTSP nodesNumber chromosomesNumber generationNumber mutationRate crossoverRate workersNumber [seed]```


## Results 
Each version of the program was tested on a Intel Xeon PHI machine from 1 to 256 thread. The tests were done with 10 generation and 20000 chromosomes and considering [500, 1000, 2000] nodes in the graph. The following speedup curves were obtained:<br/><br/>
 
<img src="https://github.com/FraCorti/TSP-GeneticAlgorithm/GeneticTSP/report/benchmark/curves/speedup_zoom_2000_20000_README.png" height="50%" width="50%"> <br/>

<img src="https://github.com/FraCorti/TSP-GeneticAlgorithm/GeneticTSP/report/benchmark/curves/speedup_zoom_1000_20000_README.png" height="50%" width="50%"><br/>

<img src="https://github.com/FraCorti/TSP-GeneticAlgorithm/GeneticTSP/report/benchmark/curves/speedup_zoom_500_20000_README.png" height="50%" width="50%">

## Acknowledgments
This project was developed for the the course of [Parallel and Distributed Systems](http://didawiki.di.unipi.it/doku.php/magistraleinformaticanetworking/spm/sdpm09support) at University of Pisa under the guide of [Prof. Marco Danelutto](http://calvados.di.unipi.it/paragroup/danelutto/) and [Prof. Massimo Torquati](http://calvados.di.unipi.it/paragroup/torquati/).