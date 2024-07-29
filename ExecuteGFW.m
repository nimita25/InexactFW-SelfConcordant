%% This code generates the input data as required and calls the function to
% implement GFW-ApproxI, GFW-ApproxII, GFW-Exact
% Change n, d, isInputDiag, epsilon, max_time parameter values
% (n,d) - Input parameters 
% isInputDiag = 1, if the input is Diag problem instance; 0 if the input
% is a Rnd problem instance
% Algorithm - Used to determine the algorithm which is executed
% Algorithm = 1 -> GFW-ApproxI, 2 -> GFW-ApproxII, 3 -> GFW-Exact
% epsilon - suboptimality error
% max_time - maximum execution time of each run of the algorithm
% number_of_runs - number of runs of the algorithms (each with a different
% seed value) so that they are randomly initialized
% initSeedValue - seed value of the first run of the algorithm. For the
% next runs, seed value is incremented by 1 each time
% The output is stored as '.mat' file in the directory 'output_GFW' within
% the current directory


% Input parameters
n = 500;
d = 100;
isInputDiag = 1;
Algorithm = 3;
epsilon = 0.05;
max_time = 72000;
number_of_runs = 10;
initSeedValue = 1;

r = n;
clear InputData;
rng(1);
if isInputDiag
    for i = 1:d
        InputData{1,i} = zeros(n,n);
        InputData{1,i}(i,i) = i;
    end
else
    for i = 1:d
        InputData{1,i} = zeros(n,n);
        for j = 1:n
            a = rand(n,1);
            InputData{1,i} = InputData{1,i} + a.*a';
        end
    end
end


GFW(InputData,isInputDiag,Algorithm,epsilon,max_time,number_of_runs,initSeedValue);
