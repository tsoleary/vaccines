%% PART 1: TOY NETWORKS
%% (1a) create toy networks: lattice, star, chain
toy_nets = ToyNets(10,10); % makes x by y lattice, N=x*y for all networks

%seed=27; %30, 32, 34, 50 for n=100 w/ 110,.025
[ER_G1]=ErdosRenyi(110,.025,randsample(1000,1));
ER = full(ER_G);
% adda 100 node ER to the toy networks
toy_nets{4} = ER;


%% (1b) run GA on each of the 4 toy networks
rng('shuffle') %random number generator seeded constantly for now

%initialize parameters
P = 50; % GA population size (change to higher val during main run)
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


% create main data structures
Toy_Solutions_Exp = {};
N_reps = 2; % number of reps
ToyMat = zeros(N*(length(transcendence_list))*length(toy_nets), 4);

counter=0; % intialize counter to index matrix

% replication loop (run for N reps and store data)
for sample=1:N_reps

    % run for a couple values of transcendence (=[2 3] shows interesting behavior)
    % store best solutions for transcendence (each net within)
    best_solutions={};

    for transcend_idx=1:length(transcendence_list)

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

        for i=1:length(toy_nets)

            %get network and N
            A=toy_nets{i};

            % run GA
            counter=counter+1; % counter to fill matrix
            
            [x,fval,exitFlag, Output] = ga(@(x) SpreadingFitnessFcnCompSize(x, A, threshold, transcendence), V, vaccineOpts);
            cur_best_solutions{i}=x;


            ToyMat(counter,1) = fval;
            ToyMat(counter,2) = Output.funccount;
            ToyMat(counter,3) = i;
            ToyMat(counter,4) = transcend_idx;

        end

        best_solutions{transcend_idx}=cur_best_solutions;
        %best_fit{transcend_idx}=cur_best_fit;
        %num_calls{transcend_idx}=cur_num_calls; % added storage of num func calls - Alex

         %look at results %(commented out while working on data structures - Alex)
        %for i=1:4
            %Visualizer(cur_best_solutions{i}, toy_nets{i}, threshold, transcendence)
        %end
    end

    Toy_Solutions_Exp{sample} = best_solutions;

end

%     Alex: (i) run this N times (20ish) to get a distribution ofsolutions/time-to-best solution (num calls).
%     (ii) Compare this to a brute-force found optimal solution (~a few hrs
%     brute force time for EACH net of size 100 choose 3

% Toy Solutions is a Matrix wit four columns (fitness, # calls, trans value, network type
% 1 2 3 and 4 (lattice, star, chain and ER respectively))

%% use best_solutions and best_fit with ER for a big figure

% Thomas plot the four best toy NWs (chain, lattice, star and er) in a 2 by
% 2 panel -Alex

% Load saved figures
s=hgload('starExample.fig');
c=hgload('chainExample.fig');
l=hgload('latticeExample.fig');
e=hgload('erExample.fig');

% Prepare subplots
figure
h(1)=subtightplot(2,2,1,[0.05,0.01], 0.075, 0.05);
set(gca,'XColor', 'none','YColor','none')
h(2)=subtightplot(2,2,2,[0.05,0.01], 0.075, 0.05);
set(gca,'XColor', 'none','YColor','none')
h(3)=subtightplot(2,2,3,[0.05,0.01]);
set(gca,'XColor', 'none','YColor','none')
h(4)=subtightplot(2,2,4,[0.05,0.01]);
set(gca,'XColor', 'none','YColor','none')

% Paste figures on the subplots
copyobj(allchild(get(s,'CurrentAxes')),h(1));
copyobj(allchild(get(c,'CurrentAxes')),h(2));
copyobj(allchild(get(l,'CurrentAxes')),h(3));
copyobj(allchild(get(e,'CurrentAxes')),h(4));

% Add legends
l(1)=title(h(1),'\fontsize{36}Star');
l(2)=title(h(2),'\fontsize{36}Chain');
l(3)=title(h(3),'\fontsize{36}Lattice');
l(4)=title(h(4),'\fontsize{36}Erd?s?R�nyi');



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

N_reps = 2; % number of replicates

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

ER_Mat = zeros(N_reps, 2); % final ER matrix

for n = 1:N_reps

    % run GA
    [x,fval, exitFlag, Output] = ga(@(x) SpreadingFitnessFcnCompSize(x, ER_G, threshold, transcendence), V, vaccineOpts);
    
    ER_Mat(n, 1) = fval;
    ER_Mat(n, 2) = Output.funccount;
    
% Visualizer(x, ER_G, threshold, transcendence)

end

% TODO: run X times to get a distribution of fitnessess (& time to best)
% completed: ER_Mat is a matrix where the first column are the fitness
% values and second are number 0f function calls. Every row is a replicate
% run

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

% completed by Thomas Above




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

for transcendence=[2] % since we don't know what it is, try a few
    for i=length(flu_nets):-1:1 %for flu networks of size ~20 up to ~1000

        flu_net=flu_nets{i};
        %N = size(flu_net,1);
        population=zeros(P,V);
        for j=1:P
            population(j,:)=randsample(1:N,V);
        end

        vaccineOpts.InitialPopulation=population;

        % RUN GA X TIMES
        for run=1:num_GA_runs

            % run GA
            [x, fval, exitFlag, Output] = ga(@(y) SpreadingFitnessFcnCompSize(y, flu_net, threshold, transcendence), V, vaccineOpts);

            Visualizer(x, flu_net, threshold, transcendence)
        end

    end
end

%TODO: Extract data for X number of runs of the above


%% (3c) Figure 3: Expected # strains in outbreak (fitness) vs network size (N), by transcendence.
% uses results of 3b



%% Example of our networks

% Load saved figures
one=hgload('real0.fig');
two=hgload('real1.fig');
three=hgload('real2.fig');
four=hgload('real3.fig');
five=hgload('real4.fig');
six=hgload('real5.fig');
seven=hgload('real7.fig');
eight=hgload('real8.fig');
nine=hgload('real9.fig');

% Prepare subplots
figure
h(1)=subtightplot(3,3,1,[0.05,0.01], 0.075, 0.05);
set(gca,'XColor', 'none','YColor','none')
h(2)=subtightplot(3,3,2,[0.05,0.01], 0.075, 0.05);
set(gca,'XColor', 'none','YColor','none')
h(3)=subtightplot(3,3,3,[0.05,0.01], 0.075, 0.05);
set(gca,'XColor', 'none','YColor','none')
h(4)=subtightplot(3,3,4,[0.05,0.01]);
set(gca,'XColor', 'none','YColor','none')
h(5)=subtightplot(3,3,5,[0.05,0.01]);
set(gca,'XColor', 'none','YColor','none')
h(6)=subtightplot(3,3,6,[0.05,0.01]);
set(gca,'XColor', 'none','YColor','none')
h(7)=subtightplot(3,3,7,[0.05,0.01]);
set(gca,'XColor', 'none','YColor','none')
h(8)=subtightplot(3,3,8,[0.05,0.01]);
set(gca,'XColor', 'none','YColor','none')
h(9)=subtightplot(3,3,9,[0.05,0.01]);
set(gca,'XColor', 'none','YColor','none')



% Paste figures on the subplots
copyobj(allchild(get(one,'CurrentAxes')),h(1));
copyobj(allchild(get(two,'CurrentAxes')),h(2));
copyobj(allchild(get(three,'CurrentAxes')),h(3));
copyobj(allchild(get(four,'CurrentAxes')),h(4));
copyobj(allchild(get(five,'CurrentAxes')),h(5));
copyobj(allchild(get(six,'CurrentAxes')),h(6));
copyobj(allchild(get(seven,'CurrentAxes')),h(7));
copyobj(allchild(get(eight,'CurrentAxes')),h(8));
copyobj(allchild(get(nine,'CurrentAxes')),h(9));

% Add legends
l(1)=title(h(1),'\fontsize{24}Star');
l(2)=title(h(2),'\fontsize{24}Chain');
l(3)=title(h(3),'\fontsize{24}Lattice');
l(4)=title(h(4),'\fontsize{24}Erd?s?R�nyi');
l(5)=title(h(5),'\fontsize{24}Star');
l(6)=title(h(6),'\fontsize{24}Chain');
l(7)=title(h(7),'\fontsize{24}Lattice');
l(8)=title(h(8),'\fontsize{24}Erd?s?R�nyi');
l(9)=title(h(9),'\fontsize{24}Erd?s?R�nyi');





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
