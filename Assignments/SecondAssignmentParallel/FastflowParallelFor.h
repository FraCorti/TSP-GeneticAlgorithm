//
// Created by checco on 23/04/20.
//

#ifndef SECONDASSIGNMENTPARALLEL__FASTFLOWPARALLELFOR_H_
#define SECONDASSIGNMENTPARALLEL__FASTFLOWPARALLELFOR_H_

#include <ff/parallel_for.hpp>
#include <mutex>

using ull = unsigned long long;

struct ParallelFor {
 private:
  const ull startPrimes;
  const ull endPrimes;
  const ull workersNumber;
  std::vector<int> primesNumber;
 public:
  ParallelFor(const ull start_primes, const ull end_primes, const ull workers_number)
      : startPrimes(start_primes), endPrimes(end_primes), workersNumber(workers_number) {
    primesNumber.reserve((size_t) (endPrimes - startPrimes) / log(startPrimes));
  }

  void Run() {
    ff::ParallelFor pf(workersNumber);
    std::mutex mutex;
    pf.parallel_for(startPrimes, endPrimes,     // start, stop indexes
                    1,                               // step size
                    0,                              // chunk size (0=static, >0=dynamic)
                    [&](const ull n) {
                      if (is_prime(n)) {
                        mutex.lock();
                        primesNumber.emplace_back(n);
                        mutex.unlock();
                      }
                    }
    );
    std::sort(primesNumber.begin(),
              primesNumber.end());      // can be parallelize with Intel Tbb installed with std::execution::par_unseq
    std::cout << "Found " << primesNumber.size() << " primes: ";
    for (auto &&prime: primesNumber) {
      std::cout << prime << ", ";
    }
    std::cout << std::endl;
  }

};
#endif //SECONDASSIGNMENTPARALLEL__FASTFLOWPARALLELFOR_H_