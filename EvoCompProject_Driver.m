% Evolutionary Computation
% Project Driver
% 10, November 2019
% Blake Williams, Thomas O'Leary and Alex Burnham

%build network
rng(9)
n=200;
A = rand(n)>.98;
A = triu(A) + triu(A,1)';
A = A - diag(diag(A));
G=graph(A);
[bin binsize]=conncomp(G);
idx = binsize(bin) == max(binsize);
GC = subgraph(G, idx);
n=GC.numnodes;
adj_mat_network=full(adjacency(GC)); %%
figure; plot(GC)
%%
% Starting parameters:
N=n;
P = 50;
global V
V = 4; % # vaccines

mutProb = 1/V; % probabilty of mutation

threshold = .1;
transcendence = 5;

% create weighted binary genome
genomeMat=zeros(P,V);
for i=1:P
    genomeMat(i,:)=randsample(1:N,V);
end

% set options for ga toolbox for bitstring
vaccineOpts = gaoptimset(...
    'PopulationType', 'doubleVector', ...
    'InitialPopulation', genomeMat, ...
    'PopulationSize',P,...
    'CrossoverFcn', @crossoversinglepoint, ...
    'CrossoverFraction', 0.2, ...
    'SelectionFcn',{@selectiontournament,2}, ...
    'Vectorized','on',...
    'Generations', 100, ...
    'FitnessLimit', 0, ...
    'MutationFcn', {@randomResetMutation, mutProb},...
    'PlotFcn',{@gaplotbestf});


% run GA (fitnessFunction, genomeLength, passOptions)
[x,fval] = ga(@(x) SpreadingFitnessFcn(x, adj_mat_network, threshold, transcendence), V, vaccineOpts);

