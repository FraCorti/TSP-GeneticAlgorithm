#!/bin/sh
clear
echo "Running benchmark "

#for ((w = 1; w <= 250; w=w+1)); do
#    echo $w;
# done
# buildFolder="cmake-build-release-local"
MW_ExecutablePath="builds/finding_primes_MW"
PF_ExecutablePath="builds/finding_primes_PF"
filenameRoot="benchmark"
maxThreads=250

for dim in 500000 5000000 50000000; do
  echo ${dim}
  for exec in ${MW_ExecutablePath} ${PF_ExecutablePath}; do
    if [ "${exec}" = "${MW_ExecutablePath}"  ]; then echo \"MasterWorker\">> ${filenameRoot}-${dim}.txt; fi
    if [ "${exec}" = "${PF_ExecutablePath}" ]; then echo \"ParallelFor\">> ${filenameRoot}-${dim}.txt; fi
    : $((w = 1))
    while [ $((w <= maxThreads)) -ne 0 ]; do
      echo $w
      ./${exec} 1 ${dim} ${w} >>${filenameRoot}-${dim}.txt
      : $((w = w + 1))
      sleep 0.3s
    done
    echo >>${filenameRoot}-${dim}.txt
    echo >>${filenameRoot}-${dim}.txt
  done
  echo
done

echo "Benchmark ended"


#for dim in 100000 1000000; do
#  echo "\"MasterWorker-${dim}\"" >>benchmarkResultsMW.txt
#  : $((w = 1))
#  while [ $((w <= maxThreads)) -ne 0 ]; do
#    echo $w
#    ./${buildFolder}/${executableName} 1 ${dim} ${w} >>benchmarkResultsMW.txt
#    : $((w = w + 1))
#    sleep 1s
#  done
#  echo >>benchmarkResultsMW.txt
#  echo >>benchmarkResultsMW.txt
#done
