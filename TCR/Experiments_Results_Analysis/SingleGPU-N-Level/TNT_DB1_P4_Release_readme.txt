data should be in the same directory as the source code, they should also have the "Result_out" folder which is the folder that holds the results for each n-nucleotide and VJ pair. The following is the required steps to compile a code and submit the job on HPC (ocelote):
1) Put all the source codes and data into one specific directory.
2) Go to that specific directory.
3) module unload gcc
4) module load cuda80
5) make 
6) modify the script base on their own directory and then
7) qsub ./Script.pbs