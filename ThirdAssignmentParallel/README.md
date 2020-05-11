# Parallel Divide and Conquer
The purpose of this assignment was to implement a parallel version of the Divide and Conquer skeleton. I started from a given sequential implementation, the template function *sequentialDc* in [this file](divideAndConquer.cpp).

## Challenges
A problem that I faced since the beginning was that a parallel version of the DaC algorithm needed to keep under control the degree of parallelism of the program. 
In fact, being DaC an algorithm that **recurses** untill the base case is reached, letting the concurrency level under no leash would mean to create hundreds, possibly thousands, of threads.  
For this reason I added a parameter called ```cutoff``` used to switch to the sequential version of the DaC algorithm from a certain point on. 
The *user* of the skeleton can then determine a good check to use as a stopping condition in the function ```checkCutoff()```.


## Example usage
To show a real usage of the algorithm I imported in the main part of the implementation of the MergeSort based on the work of [Tiziano de Matteis](https://github.com/TizianoDeMatteis).

## Building
To generate the build system do:  
```
mkdir -p build
cd build
cmake ..
```
in the root of the cloned project.  

Now build the project:
```
cmake --build .
```


## Running
The executable can be found in the ```build/``` directory. Just run it with:

```./build/ThirdAssignmentParallel``` 


## Acknowledgments
This project was developed for an assignment of the course [Parallel and Distributed Systems](http://didawiki.di.unipi.it/doku.php/magistraleinformaticanetworking/spm/sdpm09support) at University of Pisa under the guide of [Prof. Marco Danelutto](http://calvados.di.unipi.it/paragroup/danelutto/) and [Prof. Massimo Torquati](http://calvados.di.unipi.it/paragroup/torquati/).
