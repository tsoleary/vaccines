% Evolutionary Computation
% Project Driver
% 18, November 2019
% Blake Williams, Thomas O'Leary and Alex Burnham

%build network
rng(9)
n=200;
A = rand(n)>.985;
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
P = 100;
global V
V = 3; % set our threashold for what our weighted genome sums to

mutProb = 1/V; % probabilty of mutation
probZero = 0.9; % probability allele is 0

threshold = .5;
transcendence = 1;

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
    'CrossoverFraction', 0.1, ...
    'SelectionFcn',{@selectiontournament,2}, ...
    'Vectorized','on',...
    'Generations', 100, ...
    'FitnessLimit', 0, ...
    'MutationFcn', {@randomResetMutationFcn, mutProb},...
    'PlotFcn',{@gaplotbestf});


% run GA (fitnessFunction, genomeLength, passOptions)
[x,fval] = ga(@(x) SpreadingFitnessFcn(x, adj_mat_network, threshold, transcendence), V, vaccineOpts);

