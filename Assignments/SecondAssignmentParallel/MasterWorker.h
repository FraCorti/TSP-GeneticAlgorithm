//
// Created by checco on 23/04/20.
//

#ifndef SECONDASSIGNMENTPARALLEL__MASTERWORKER_H_
#define SECONDASSIGNMENTPARALLEL__MASTERWORKER_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ff/parallel_for.hpp>
#include <ff/pipeline.hpp>
#include <ff/farm.hpp>
#include <tuple>

using ull = unsigned long long;

// see http://en.wikipedia.org/wiki/Primality_test
static bool is_prime(ull n) {
  if (n <= 3) return n > 1; // 1 is not prime !
  if (n % 2 == 0 || n % 3 == 0) return false;
  for (ull i = 5; i * i <= n; i += 6) {
    if (n % i == 0 || n % (i + 2) == 0)
      return false;
  }
  return true;
}

//! worker node
struct Worker : ff::ff_node_t<std::pair<ull, ull>, std::vector<ull>> {

  std::vector<ull> *svc(std::pair<ull, ull> *data) override {
    ull indexStart = data->first;
    ull indexEnd = data->second;
    delete data;
    auto *primesNumber = new std::vector<ull>;
    primesNumber->reserve((size_t) (indexEnd - indexStart) / log(indexStart));

    // iterating over the range and test primality
    for (; indexStart <= indexEnd; indexStart++) {
      if (is_prime((indexStart)))
        primesNumber->emplace_back(indexStart);
    }
    return primesNumber;
  }
};

//! master node
struct Master : ff::ff_node_t<std::vector<ull>, std::tuple<ull>> {
 private:
  const ull startPrimes;
  const ull endPrimes;
  const ull workersNumber;
  size_t nWorkers = 0;
  std::vector<ull> primesFound;
 public:
  Master(ull start_primes,
         ull end_primes,
         ull workers_Number)
      : startPrimes(start_primes), endPrimes(end_primes), workersNumber(workers_Number) {
    primesFound.reserve((size_t) (endPrimes - startPrimes) / log(startPrimes));
  }

  std::tuple<ull> *svc(std::vector<ull> *primeNumbers) override {

    // divide works among workers and push to the code the tuple with primes ranges
    if (primeNumbers == nullptr) {
      ull workerCellsToCompute =
          std::floor(static_cast<float>(endPrimes - startPrimes) / static_cast<float>(workersNumber));
      ull remainedElements = (endPrimes - startPrimes) % workersNumber;
      ull indexStart = startPrimes;
      ull indexEnd;

      for (size_t iteration = 0; iteration < workersNumber; iteration++) {
        indexEnd = indexStart + workerCellsToCompute;
        ff_send_out(new std::pair(indexStart, indexEnd));
        indexStart = indexEnd++;
      }

      // consider extra elements
      if (remainedElements) {
        indexEnd++;
        for (; indexEnd <= endPrimes; indexEnd++) {
          if (is_prime(indexEnd))
            primesFound.push_back(indexEnd);
        }
      }
      return GO_ON;
    }

    // manage returned primes from workers
    std::vector<ull> &result = *primeNumbers;
    primesFound.insert(primesFound.end(), result.begin(), result.end());
    delete primeNumbers;
    if (++nWorkers == workersNumber) return EOS;
    return GO_ON;
  }

  void svc_end() override {
    std::sort(primesFound.begin(),
              primesFound.end());
    std::cout << "Found " << primesFound.size() << " primes: ";
    for (auto &&prime: primesFound) {
      std::cout << prime << ", ";
    }
    std::cout << std::endl;
  }
};
#endif //SECONDASSIGNMENTPARALLEL__MASTERWORKER_H_