function [] = GFW(InputData, isDiagInput, AlgorithmType, epsilon, max_time, number_of_runs,initSeedValue)

%% Algorithm to solve min - \sum_{i=1}^d log(<A_i,X>) s.t. Tr(X) = 1, X p.s.d.
% The algorithm does not store the decision variable X, it tracks the
% progress of the algorithm via the d-dimensional vector vector v, where
% v_i = <A_i,X> for all i = 1,...,d
% The algorithm accepts two types of input instances: Diagonal, Random
% Diagonal: For each i, i-th diagonal entry of A_i = i, with all other
% entries 0
% Random: For each i, A_i is a sum of j rank-1 matrices (u_j)(u_j)^T such
% that each entry of (u_j) ~ N(0,1)
% Input arguments:
% InputData - A sctructure array with d entries such that each entry is an
% nxn matrix A_i
% isDiagInput - takes value 0 or 1. If 1, A_i's are diagonal matrices as
% described in the paper; if 0, A_i's are random matrices
% AlgorithmType - takes values 1,2,3 for GFW_ApproxI, GFW_ApproxII, and
% GFW_Exact algorithms respectively
% epsilon - suboptimality error (stopping criteria) of the algorithm
% max_time - maximum execution time of each run of the algorithm
% number_of_runs - number of runs of the algorithms (each with a different
% seed value) so that they are randomly initialized
% initSeedValue - seed value of the first run of the algorithm. For the
% next runs, seed value is incremented by 1 each time
% The output is stored as '.mat' file in the directory 'output_GFW' within
% the current directory


%% Check input data
narginchk(6,7);
if nargin == 6
    max_time = 36000;
end

%% Load input data
A = InputData;
clear InputData;

%% Set parameters
%nsamples = 1; %number of Gaussian random vectors z ~ N(0,X) stored
d = length(A);
n = length(A{1,1});
if AlgorithmType == 3
    etak = 0;
else
    etak = epsilon/2;
end


%% Start the runs
for nr = 1:number_of_runs
    %% Set seed value
    rng(initSeedValue+nr-1);

    disp(['Starting run ',int2str(nr)]);

    %% Initialize the variables
    disp('Initializing the algorithm');
    v = zeros(d,1);
    traceX = 0;
    for ii = 1:d
        Xfact = randn(n,1);
        for kk = 1:d
            v(kk) = v(kk)+Xfact'*A{1,kk}*Xfact;
        end
        traceX = traceX + norm(Xfact,2)^2;
    end
    v = v/traceX; % Define v_i = Tr(A_i^TX) for all i
    InitialFnVal = -sum(log(v));


    %% Compute instance-specific parameters for given initialization
    % Compute delta0 (initial optimality gap)
    if isDiagInput % When [A_i]_{ii} = i
        tv = 1:d;
        optVal = -sum(log(tv/d));
        delta0 = -sum(log(v))+sum(log(tv/d));
    else
        delta0 = d*log(n); % For random A_i matrices
    end

    % Compute upper bounds on KC1, KC2, KC2Gk
    % The upper bounds are KC1u, KC2u, KC2Gku which are defined as
    % K_{C1}^u, K_{Delta}^u, K_{G_t^a}^u in the paper
    if AlgorithmType == 1
        KC1u = ceil(10.6*delta0);
        KC2u = ceil(12*d^2*max(0,((1/(epsilon-etak))-(1/(delta0-etak)))));
    elseif AlgorithmType == 2
        KC1u = ceil(10.6*delta0);
        KC2u = ceil(192*d^2*((1/epsilon)-(1/(delta0))))+ ceil(log2((5*d)/(4*epsilon))+1);
    else
        KC1u = ceil(5.3*(delta0+d)*log(10.6*delta0));
        KC2u = ceil(12*d^2*max(0,((1/epsilon)-(1/delta0))));
    end
    KC2Gku = ceil(24*d^2/epsilon);

    %% Approximate Generalized Frank-Wolfe
    tic;
    disp('Starting the algorithm...');

    %% First iteration of the algorithm
    %Compute the gradient of obj function
    gradf = sparse(n,n);
    for ii = 1:d
        gradf = gradf - A{1,ii}/v(ii);
    end

    %Generate the update direction using 'eigs'
    %opts.maxit = 400;
    if AlgorithmType == 3
        [U,L] = eig(gradf);
        [l,pos] = min(diag(L));
        u = U(:,pos);
    elseif AlgorithmType == 2
        norm_gradf = norm(gradf,2);
        opts.tol = min(etak/norm_gradf,1);
        [u,l] = eigs(gradf,1,'SR',opts);
    else
        [u,l] = eigs(gradf,1,'SR');
    end

    %Compute update direction
    h = zeros(d,1);
    if l <= 0
        for ii = 1:d
            h(ii) = u'*A{1,ii}*u;
        end
    end

    % Compute Frank-Wolfe duality gap
    Gk = (-1./v)'*(v-h);
    Dk = sqrt(((h-v)./(v.^2))'*(h-v));
    gamma = (Gk)/(Dk*((Gk)+Dk));
    gammak = min(1,gamma);
    if isDiagInput
        deltak = -sum(log(v)) - optVal;
    end
    minGk = Gk;
    
    t = 1; %Set iteration counter
    % Initialize arrays to track the change in G_k^a, function value
    GktrackC1 = [];
    GktrackC2 = [];
    FnValTrack = [];
    
    % Initialize the counter of number of iterations of each case
    %C1 - G_k^a > d, C2 - G_k^a <= d & delta_k (optimality gap) >= epsilon
    %C3 - epsilon <= G_k^a <= d
    NItrC1 = 0;
    NItrC2 = 0;
    NItrC3 = 0;
    
    % Track change in optimality gap for Diag problem instances
    if isDiagInput
        deltatrackC1 = [];
        deltatrackC2 = [];
    end

    %% Run the loop until deltak (optimality gap) <= epsilon
    if isDiagInput
        disp('Reducing optimality gap to <= epsilon');
        disp('#itr|Duality Gap');
        if Gk > d
            GktrackC1 = [GktrackC1,Gk];
            deltatrackC1 = [deltatrackC1,deltak];
        else
            GktrackC2 = [GktrackC2,Gk];
            deltatrackC2 = [deltatrackC2,deltak];
        end

        while ( deltak > epsilon && toc <= max_time )
            %% Display intermediate status
            if mod(t,100) == 0
                disp([int2str(t),'|',num2str(round(Gk,2))]);
            end

            %% Update the linear map 'v'
            v = (1-gammak)*v + gammak*h;
            t = t+1; % Increment iteration counter

            %% Compute the update direction using 'eigs'
            %Compute gradient of obj function at the new point
            gradf = zeros(n,n);
            for ii = 1:d
                gradf = gradf - A{1,ii}/v(ii);
            end
            %Generate the update direction using 'eigs'
            if AlgorithmType == 3
                [U,L] = eig(gradf);
                [l,pos] = min(diag(L));
                u = U(:,pos);
            elseif AlgorithmType == 2
                norm_gradf = norm(gradf,2);
                etak = (1/4)*minGk;
                opts.tol = min(etak/norm_gradf,1);
                [u,l] = eigs(gradf,1,'SR',opts);
            else
                [u,l] = eigs(gradf,1,'SR');
            end

            %Compute update direction
            h = zeros(d,1);
            if l <= 0
                for ii = 1:d
                    h(ii) = u'*A{1,ii}*u;
                end
            end

            % Compute Frank-Wolfe duality gap
            Gk = (-1./v)'*(v-h);
            deltak = -sum(log(v)) - optVal;
            minGk = min(Gk,minGk);
            Dk = sqrt(((h-v)./(v.^2))'*(h-v));
            gamma = (Gk)/(Dk*((Gk)+Dk));
            gammak = min(1,gamma); % Determine step length
            
            % Update the arrays to track the duality gap, optimality gap
            % Increment the iteration counter for each case
            if Gk > d
                GktrackC1 = [GktrackC1,Gk];
                deltatrackC1 = [deltatrackC1,deltak];
                NItrC1 = NItrC1 + 1;
            else
                GktrackC2 = [GktrackC2,Gk];
                deltatrackC2 = [deltatrackC2,deltak];
                NItrC2 = NItrC2 + 1;
                NItrC3 = NItrC3 + 1;
            end
        end
        
        if toc <= max_time
            disp('Epsilon optimal solution found. Optimality gap <= epsilon');
            %disp(-sum(log(v)));
        end
    end

    %% Run the loop until Gk <= epsilon
    if toc <= max_time
        disp('Reducing duality gap to <= epsilon');
        disp('#itr|Duality Gap');

        while ( Gk > epsilon && toc <= max_time )
            %% Display intermediate status
            if mod(t,100) == 0
                disp([int2str(t),'|',num2str(round(Gk,2))]);
            end

            %% Update the linear map 'v'
            v = (1-gammak)*v + gammak*h;
            t = t+1; % Increment iteration counter
            
            %% Compute the update direction using 'eigs'
            %Compute gradient of obj function at the new point
            gradf = zeros(n,n);
            for ii = 1:d
                gradf = gradf - A{1,ii}/v(ii);
            end
            %Generate the update direction using 'eigs'
            if AlgorithmType == 3
                [U,L] = eig(gradf);
                [l,pos] = min(diag(L));
                u = U(:,pos);
            elseif AlgorithmType == 2
                norm_gradf = norm(gradf,2);
                etak = (1/4)*minGk;
                opts.tol = min(etak/norm_gradf,1);
                [u,l] = eigs(gradf,1,'SR',opts);
            else
                [u,l] = eigs(gradf,1,'SR');
            end
            
            %Compute update direction
            h = zeros(d,1);
            if l <= 0
                for ii = 1:d
                    h(ii) = u'*A{1,ii}*u;
                end
            end

            % Compute Frank-Wolfe duality gap
            Gk = (-1./v)'*(v-h);
            minGk = min(Gk,minGk);
            Dk = sqrt(((h-v)./(v.^2))'*(h-v));
            gamma = (Gk)/(Dk*((Gk)+Dk));
            gammak = min(1,gamma); % Determine step length

            
            % Update the arrays to track the duality gap
            % Increment the iteration counter for each case
            if Gk > d
                GktrackC1 = [GktrackC1,Gk];
                NItrC1 = NItrC1+1;
            else
                GktrackC2 = [GktrackC2,Gk];
                NItrC3 = NItrC3 + 1;
            end
            FnValTrack = [FnValTrack,-sum(log(v))];
        end
    end

    %% Generate objective function value
    fnval = -sum(log(v));
    time = toc;

    %% Display output
    if time > max_time
        disp('Maximum time reached. Solution not found');
    elseif Gk <= epsilon
        disp('Epsilon-optimal solution found');
    end
    disp('n,d');
    disp([n,d]);
    disp('Number of iterations of Case I (Gk > d)');
    disp(NItrC1);
    if isDiagInput
        disp('Number of iterations of Case II (deltak > epsilon & Gk <= d)');
        disp(NItrC2);
    end
    disp('Number of iterations of Case III (Gk <= epsilon)');
    disp(NItrC3);
    disp('Time required');
    disp(time);
    disp('Duality gap');
    disp(Gk);
    disp('Objective function value');
    disp(fnval);
    %disp('Gaussian samples with covariance equal to the solution generated. You can use them to recreate the p.s.d. d.v. using covariance estimation techniques');

    %% Write output
    % Input parameters
    SCLH.InputParams.n = n;
    SCLH.InputParams.d = d;
    SCLH.InputParams.epsilon = epsilon;
    if isDiagInput
        SCLH.InputParams.optVal = optVal;
    end
    SCLH.InputParams.delta0{nr} = delta0;
    SCLH.InputParams.InitialFnVal{nr} = InitialFnVal;
    SCLH.InputParams.KUpperBoundC1{nr} = KC1u;
    SCLH.InputParams.KUpperBoundC2{nr} = KC2u;
    SCLH.InputParams.KGkConvUpperBoundC2 = KC2Gku;
    if AlgorithmType == 1
        SCLH.InputParams.Algorithm = 'Approximate algorithm with constant eta';
        Algo = 'ApproxI';
    elseif AlgorithmType == 2
        SCLH.InputParams.Algorithm = 'Approximate algorithm with varying eta';
        Algo = 'ApproxII';
    else
        SCLH.InputParams.Algorithm = 'Exact algorithm';
        Algo = 'Exact';
    end
    RunNo = ['Run',int2str(nr)];
    SeedValue = initSeedValue+nr-1;
    SCLH.InputParams.Initialization{nr} = ['Random initialization with N(0,1) for each entry and seed value ',int2str(SeedValue)];
    
    % Output parameters
    SCLH.OutputParams.KC1{nr} = NItrC1;
    if isDiagInput
        SCLH.OutputParams.KC2{nr} = NItrC2;
    else
        SCLH.OutputParams.UpdatedDelta0{nr} = InitialFnVal - fnval + epsilon;
        delta0 = InitialFnVal - fnval + epsilon;
        if AlgorithmType == 1
            KC1u = ceil(10.6*delta0);
            KC2u = ceil(12*d^2*max(0,((1/(2*epsilon-etak))-(1/(delta0-etak)))));
        elseif AlgorithmType == 2
            KC1u = ceil(10.6*delta0);
            KC2u = ceil(192*d^2*((1/(2*epsilon))-(1/(delta0))))+ ceil(log2((5*d)/(8*epsilon))+1);
        else
            KC1u = ceil(5.3*(delta0+d)*log(10.6*delta0));
            KC2u = ceil(12*d^2*max(0,((1/(2*epsilon))-(1/delta0))));
        end
        SCLH.OutputParams.UpdatedKC1u{nr} = KC1u;
        SCLH.OutputParams.UpdatedKC2u{nr} = KC2u;
        FnValTrack = sort(FnValTrack,'descend');
        FnValTrack = FnValTrack<=fnval+epsilon;
        SCLH.OutputParams.ApproxKC2{nr} = find(FnValTrack,1,'first')-NItrC1;
    end
    SCLH.OutputParams.KGkConvC2{nr} = NItrC3;
    SCLH.OutputParams.TotalTimeSecs{nr} = time;
    SCLH.OutputParams.FinalObjFnVal{nr} = fnval;
    SCLH.OutputParams.Gk{nr} = Gk;
    if isDiagInput
        SCLH.OutputParams.deltatrackC1{nr} = deltatrackC1;
        SCLH.OutputParams.deltatrackC2{nr} = deltatrackC2;
    end
end

if ~exist('output_GFW','dir'), mkdir('output_GFW'); end
filename = [datestr(now,'dd-mm-yy-HH:MM-'),Algo,'-',int2str(n),'-',int2str(d)];
%filename = [datestr(now,'dd-mm-yy-HH:MM-'),Algo,'R1','-',int2str(n),'-',int2str(d)];
save(['output_GFW/',filename],'SCLH','-v7.3');
