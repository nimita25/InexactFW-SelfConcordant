ExecuteGFW.m generates the input data as required and calls the function to implement GFW-ApproxI, GFW-ApproxII, GFW-Exact

In this file (ExecuteGFW.m), we
1) Define the input parameters,
2) Generate the input problem instance according to the input data size (n, d)
3) Call the function GFW which executes the algorithm.

The function GFW is defined in the file GFW.m. Please ensure that ExecuteGFW.m and GFW.m are stored in the same directory.

To execute any algorithm, open the file ExecuteGFW.m and change the value of the following parameters (on lines 20-27):
1) (n,d) - Input parameters that determine the size of the input problem instance
2) isInputDiag = 1, if the input is Diag problem instance; 0 if the input is a Rnd problem instance (please refer to the paper for more details about the instances)
3) Algorithm - Used to determine the algorithm which is executed. Algorithm = 1 -> GFW-ApproxI, 2 -> GFW-ApproxII, 3 -> GFW-Exact
4) epsilon - suboptimality error of the algorithm
5) max_time - maximum execution time of each run of the algorithm
6) number_of_runs - number of runs of the algorithms (each with a different seed value) so that they are randomly initialized
7) initSeedValue - seed value of the first run of the algorithm. For the next runs, seed value is incremented by 1 each time


The output is stored as '.mat' file in the directory 'output_GFW' within the current directory
