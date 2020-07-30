#!/bin/sh
clear
echo "Running benchmark "

FF_ExecutablePath="builds/TSP_Fastflow"
Standard_ExecutablePath="builds/TSP_Standard"
filenameRoot="benchmark"
maxThreads=250
mutationRate=0.1
crossoverRate=0.1
generationsNumber=10

for nodes in 500 1000 2000; do
  echo ${nodes}
  for chromosomes in 500 5000 20000; do
    for exec in ${FF_ExecutablePath} ${Standard_ExecutablePath}; do
      if [ "${exec}" = "${FF_ExecutablePath}"  ]; then echo \"Fastflow\">> ${filenameRoot}-${nodes}-${chromosomes}.txt; fi
      if [ "${exec}" = "${Standard_ExecutablePath}" ]; then echo \"Standard\">> ${filenameRoot}-${nodes}-${chromosomes}.txt; fi
      : $((w = 1))
      while [ $((w <= maxThreads)) -ne 0 ]; do
        echo $w
        ./${exec} ${nodes} ${chromosomes} ${generationsNumber} ${mutationRate} ${crossoverRate} ${w} >>${filenameRoot}-${nodes}-${chromosomes}.txt
        : $((w = w + 1))
      done
      echo >>${filenameRoot}-${nodes}-${chromosomes}.txt
      echo >>${filenameRoot}-${nodes}-${chromosomes}.txt
    done
  done
  echo
done
echo "Benchmark ended"