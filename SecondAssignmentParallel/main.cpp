#include "MasterWorker.h"
#include "FastflowParallelFor.h"

int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cerr << "use: " << argv[0] << " number1 number2 workersNumber chunkSize [print=off|on]\n";
    return -1;
  }
  const ull n1 = std::stoll(argv[1]);
  const ull n2 = std::stoll(argv[2]);
  const ull workersNumber = std::stoll(argv[3]);

  // create master
  Master master(n1, n2, workersNumber);

  // create workers
  std::vector<std::unique_ptr<ff::ff_node>> workers;
  for (ull i = 0; i < workersNumber; i++) {
    workers.push_back(std::make_unique<Worker>());
  }

  ff::ff_Farm<> farm(std::move(workers), master);
  farm.remove_collector();
  farm.wrap_around();
  // farm.set_scheduling_ondemand();

  ff::ffTime(ff::START_TIME);
  if (farm.run_and_wait_end() < 0) {
    ff::error("running farm");
    return -1;
  }
  ff::ffTime(ff::STOP_TIME);
  std::cout << "Master workers time: " << ff::ffTime(ff::GET_TIME) << " (ms)\n";

  ff::ffTime(ff::START_TIME);
  ParallelFor parallelFor(n1, n2, workersNumber);
  parallelFor.Run();
  ff::ffTime(ff::STOP_TIME);
  std::cout << "Parallel for time: " << ff::ffTime(ff::GET_TIME) << " (ms)\n";
  return 0;
}