#include <vector>
#include <thread>
#include <functional>
#include <iostream>

template<typename Tin, typename Tout>
Tout sequentialDc(Tin input,
           bool (*basecase)(Tin),
           Tout (*solve)(Tin),
           std::vector<Tin> (*divide)(Tin),
           Tout (*conquer)(std::vector<Tout>)) {

  if (basecase(input)) {
    return (solve(input));
  } else {
    auto subproblems = divide(input);

    std::transform(subproblems.begin(),
                   subproblems.end(),
                   subproblems.begin(),
                   [&](Tin x) {
                     auto res = sequentialDc(x, basecase, solve, divide, conquer);
                     return (res);
                   });

    auto result = conquer(subproblems);
    return (result);
  }
};

template<typename Tin, typename Tout>
Tout parallelDc(Tin input,
                bool (*basecase)(Tin),
                Tout (*solve)(Tin),
                std::vector<Tin> (*divide)(Tin),
                Tout (*conquer)(std::vector<Tout>),
                bool (*checkCutOff)(Tin)) {
  if (basecase(input)) {
    return (solve(input));
  } else {
    // divide the input and store in subproblems
    std::vector<Tin> subproblems = divide(input);

    // initialize the vector to store the result
    std::vector<Tout> results(subproblems.size());
    int subProblemsNumber = subproblems.size();

    // list of threads
    std::vector<std::thread> threadList;

    // check granularity of the subproblems
    if(!checkCutOff(subproblems.front())){

      // declare a lambda that recurse on parallelDc
      auto threadLambda = [&results, basecase, solve, divide, conquer, checkCutOff, &subproblems](const int i) {
        results[i] = parallelDc(subproblems[i], basecase, solve, divide, conquer, checkCutOff);
      };

      for (int i = 0; i != subProblemsNumber; i++) {
        threadList.emplace_back(std::thread(threadLambda, i));
      }
    }else {

      // declare a lambda that recurse on sequentialDc
      auto threadLambda = [&results, basecase, solve, divide, conquer, checkCutOff, &subproblems](const int i) {
        results[i] = sequentialDc(subproblems[i], basecase, solve, divide, conquer);
      };

      for (int i = 0; i != subProblemsNumber; i++) {
        threadList.emplace_back(std::thread(threadLambda, i));
      }
    }

    // wait the threads
    std::for_each(threadList.begin(), threadList.end(),
                  [](std::thread & currentThread)
                  {
                    if (currentThread.joinable())
                    {
                      currentThread.join();
                    }
                  }
    );

    auto result = conquer(subproblems);
    return (result);
  }
};