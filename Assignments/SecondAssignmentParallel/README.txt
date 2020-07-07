Il programma presenta un main file (main.cpp) dove vengono incluse le classi MasterWorker.h e FastFlowParallelFor.h . Compilando il main tramite g++ e passando i parametri IndexStart IndexEnd WorkerNumbersr il programma stampa in output la quantità di numeri primi trovati, i numeri primi trovati in ordine crescente e il tempo, in millisecondi, ottenuto utilizzando Master-worker e ParallleFor.  

All'interno della cartella benchmark sono presenti i dati ottenuti sullo Xeon-Phi, òl programma è stato testato con range di 500k, 5mln e 10mln in entrambe le configurazioni (ParallelFor e MasterWorker). Per queste tre configurazioni sono presenti i grafici di speedup e scalability al variare del numero di worker all'interno della cartella /benchmark/img . Per ogni grafico è presente una versione zommata chiamata *_zoom.png . 

Il codice prodotto è disponibile su GitHub al seguente link: https://github.com/FraCorti/Parallel-and-Distributed-Systems/tree/master/SecondAssignmentParallel . È inoltre presente un README dove viene illustrato come eseguire il programma tramite CMake ed ottenere i dati tramite l'utilizzo dello script benchmark.sh
 
