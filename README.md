ExecuteGFW.m generates the input data as required and calls the function to implement GFW-ApproxI, GFW-ApproxII, GFW-Exact

In the code, change the following: n, d, isInputDiag, Algorithm, epsilon, max_time, number_of_runs, initSeedValue
1. (n,d) - Input parameters 
2. isInputDiag = 1, if the input is Diag problem instance; 0 if the input is a Rnd problem instance
3. Algorithm - Used to determine the algorithm which is executed. Algorithm = 1 -> GFW-ApproxI, 2 -> GFW-ApproxII, 3 -> GFW-Exact
4. epsilon - suboptimality error
5. max_time - maximum execution time of each run of the algorithm
6. number_of_runs - number of runs of the algorithms (each with a different seed value) so that they are randomly initialized
7. initSeedValue - seed value of the first run of the algorithm. For the next runs, seed value is incremented by 1 each time

The output is stored as '.mat' file in the directory 'output_GFW' within the current directory
