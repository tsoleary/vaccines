%% PART 1: TOY NETWORKS 
%% (1a) create toy networks: lattice, star, chain
toy_nets = ToyNets(10,10); % makes x by y lattice, N=x*y for all networks



%% (1b) run GA on each of the toy networks
rng('shuffle') %random number generator seeded constantly for now

%initialize parameters
P = 200; % GA population size
N = size(toy_nets{1},1); % # nodes
nGen=40; % GA generations, very small # needed in part 1
global V 
V = 3; % # vaccines
mutProb = 1/V; % probabilty of mutation
if V==1 % in case of 1 vaccine, set mutation=1/2 (not 1/1)
    mutProb=0.5; 
end
threshold = .5; % from 0 to 1, kept 1/2.
transcendence_list = [2 3];

%initialize population
population=zeros(P,V);
for j=1:P
    population(j,:)=randsample(N,V);
end

% run for a couple values of transcendence (=[2 3] shows interesting behavior)
% store best solutions for transcendence (each net within)
best_solutions={};
best_fit={};
for transcend_idx=length(transcendence_list)
    
    %get current transcendence value
    transcendence=transcendence_list(transcend_idx);
    % set GA options
    vaccineOpts = gaoptimset(...
        'PopulationType', 'doubleVector', ...
        'InitialPopulation', population, ...
        'PopulationSize',P,...
        'CrossoverFcn', @crossoversinglepoint, ...
        'CrossoverFraction', 0.5, ...
        'SelectionFcn',{@selectiontournament,2}, ...
        'Vectorized','on',...
        'Generations', nGen, ...
        'FitnessLimit', 0, ...
        'PlotFcn',{@gaplotbestf},...
        'MutationFcn',{@randomResetMutation, mutProb});
    % loop through each toy network
    % store best solution for each toy network
    cur_best_solutions={};
    cur_best_fit={};
    for i=1:length(toy_nets)

        %get network and N
        A=toy_nets{i};
        
        % run GA 
        [x,fval] = ga(@(x) SpreadingFitnessFcnCompSize(x, A, threshold, transcendence), V, vaccineOpts);
        cur_best_solutions{i}=x;
        cur_best_fit{i}=fval;
    end
    
    best_solutions{transcend_idx}=cur_best_solutions;
    best_fit{transcend_idx}=cur_best_fit;
    
    %  look at results
    for i=1:3
        Visualizer(cur_best_solutions{i}, toy_nets{i}, threshold, transcendence)
    end
end

%TODO (i) run this X times (20ish) to get a distribution ofsolutions/time-to-best solution. 
%     (ii) Compare this to a brute-force found optimal solution (~a few hrs
%     brute force time for EACH net of size 100 choose 3
%TODO use best_solutions and best_fit in part 2d with ER for a big figure

%% (1c) Figure 1: Diagram of fitness calculation
%TODO
%probably outside of Matlab, Blake will do this


%NOTE: if we are to brute force search any of th

%% PART 2: ERDOS RENYI RANDOM GRAPH
%% (2a) create ER random graphs
%seed=27; %30, 32, 34, 50 for n=100 w/ 110,.025
%[ER_G, n, m]=ErdosRenyi(110, .02, seed);
%[ER_G, n, m]=ErdosRenyi(0,.012,1);
[ER_G, n, m]=ErdosRenyi(150,.012,randsample(1000,1));
%% (2b) run GA on the ER graph
rng('shuffle')
%initialize
P = 200; % GA population size
N = size(ER_G,1);
nGen=30;
global V
V = 4; % # vaccines
mutProb = 1/V; % probabilty of mutation
if V==1
    mutProb=0.5;
end
threshold = .5;
transcendence = 2;
population=zeros(P,V);
for j=1:P
    population(j,:)=randsample(1:N,V);
end

% set options for ga toolbox for bitstring
vaccineOpts = gaoptimset(...
    'PopulationType', 'doubleVector', ...
    'InitialPopulation', population, ...
    'PopulationSize',P,...
    'CrossoverFcn', @crossoversinglepoint, ...
    'CrossoverFraction', 0.5, ...
    'SelectionFcn',{@selectiontournament,2}, ...%{@selectionstochunif}, ..
    'Vectorized','on',...
    'Generations', nGen, ...
    'FitnessLimit', 0, ...
    'MutationFcn', {@randomResetMutation, mutProb},...
    'PlotFcn',{@gaplotbestf});

% run GA 
[x,fval] = ga(@(x) SpreadingFitnessFcnCompSize(x, ER_G, threshold, transcendence), V, vaccineOpts);

Visualizer(x, ER_G, threshold, transcendence)

% TODO: run X times to get a distribution of fitnessess (& time to best)

%% (2c) Test random vaccination strategies
%TODO: clean this section up to find the expected best fitness for a given 
% number of random vaccinations (ie a control, to see how good the GA is).
t=6;
random_vaccs=zeros(P*t,V);
for j=1:P*t
    random_vaccs(j,:)=randsample(1:N,V);
end
trials=SpreadingFitnessFcnCompSize(random_vaccs, ER_G, threshold, transcendence);

figure;hist(trials)
title([num2str(mean(trials)) '  ' num2str(min(trials))])


%% (2d) Figure 2: Vaccinations in toy networks (verification/explanatory) 
% Show a run and explain simple happenings

figure;
% (a) Chain
subplot(1,4,1)


% (b) Star
subplot(1,4,2)


% (c) Lattice
subplot(1,4,3)


% (d) Erdos Renyi
subplot(1,4,4)






%% PART 3: RUN ON FLU NETS OF VARIOUS SIZES FOR ~3 DIFFERENT TRANSCENDENCE VALUES
%% (3a) import flu networks
load('H3N2_flu_net_components.mat'); 
flu_nets=adjmats; % ordered from largest to smallest
clear adjmats



%% (3b) run GA on flu nets

num_GA_runs=1;
P = 50; % GA population size
nGen=20;
global V
V = 4; % # vaccines
mutProb = 1/V; % probabilty of mutation
if V==1
mutProb=0.5;
end
threshold = .5;


% set options for ga toolbox for bitstring
vaccineOpts = gaoptimset(...
'PopulationType', 'doubleVector', ...
'PopulationSize',P,...
'CrossoverFcn', @crossoversinglepoint, ...
'CrossoverFraction', 0.5, ...
'SelectionFcn',{@selectiontournament,2}, ... {@selectionstochunif}
'Vectorized','on',...
'Generations', nGen, ...
'FitnessLimit', 0, ...
'MutationFcn', {@randomResetMutation, mutProb},...
'PlotFcn',{@gaplotbestf});

for transcendence=[1 1.5 2] % since we don't know what it is, try a few
    for i=length(flu_nets):-1:1 %for flu networks of size ~20 up to ~1000
        
        flu_net=flu_nets{i};
        N = size(flu_net,1);
        population=zeros(P,V);
        for j=1:P
            population(j,:)=randsample(1:N,V);
        end
        
        vaccineOpts.InitialPopulation=population;
        
        % RUN GA X TIMES
        for run=1:num_GA_runs

            % run GA 
            [x,fval] = ga(@(y) SpreadingFitnessFcnCompSize(y, flu_net, threshold, transcendence), V, vaccineOpts);

            Visualizer(x, flu_net, threshold, transcendence)
        end
            
    end
end

%TODO: Extract data for X number of runs of the above


%% (3c) Figure 3: Expected # strains in outbreak (fitness) vs network size (N), by transcendence.
% uses results of 3b






%%  PART 4: RUN ON GROWING NETWORK FLU
%% (4a) identify test network and optimal solution
%TODO: this.
%TODO: Blake provides network
%TODO



%% (4b) Grow network and evaluate fitness at each step
%TODO
%TODO
%TODO



%% (4c) Figure 4: Expected # strains in outbreak (fitness) w/network growth (N)
%TODO used data from 4b
%TODO
%TODO



%% (4d) repeat with simulated network via Alex's lego algorithm
%TODO
%TODO
%TODO