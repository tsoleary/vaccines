%% PART 1: TOY NETWORKS
%% (1a) create toy networks: lattice, star, chain
toy_nets = ToyNets(10,10); % makes x by y lattice, N=x*y for all networks

%seed=27; %30, 32, 34, 50 for n=100 w/ 110,.025
[ER_G1]=ErdosRenyi(110,.025,randsample(1000,1));
ER = full(ER_G1);
% adda 100 node ER to the toy networks
toy_nets{4} = ER;


%% (1b) run GA on each of the 4 toy networks
tic
rng('shuffle') %random number generator seeded constantly for now

%initialize parameters 100 abd 12 orignal run
P = 300; % GA population size (change to higher val during main run)
N = size(toy_nets{1},1); % # nodes
nGen=50; % GA generations, very small # needed in part 1
global V
V = 3; % # vaccines
disp(V)
mutProb = 1/V; % probabilty of mutation
if V==1 % in case of 1 vaccine, set mutation=1/2 (not 1/1)
    mutProb=0.5;
end
threshold = .5; % from 0 to 1, kept 1/2.
transcendence_list = [2 2.5 3];

%initialize population
population=zeros(P,V);
for j=1:P
    population(j,:)=randsample(N,V);
end


% create main data structures
Toy_Solutions_Exp = {};
N_reps = 20; % number of reps
ToyMat = zeros(N_reps*(length(transcendence_list))*length(toy_nets), 4);

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
            'PlotFcn',{@gaplotbestf },...''
            'MutationFcn',{@randomResetMutation, mutProb},....
            'OutputFcn',{@outputFcn});

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
            record = outputFcn();
            fcncalls_to_best_solution = P * (record(end).LastImprovement + 1)
            clear outputFcn

            ToyMat(counter,1) = fval;
            %ToyMat(counter,2) = Output.funccount;
            %changed to: time to best solution
            ToyMat(counter,2) = fcncalls_to_best_solution;
            ToyMat(counter,3) = i;
            ToyMat(counter,4) = transcendence;

        end

        best_solutions{transcend_idx}=cur_best_solutions;
        %best_fit{transcend_idx}=cur_best_fit;
        %num_calls{transcend_idx}=cur_num_calls; % added storage of num func calls - Alex

         %look at results %(commented out while working on data structures - Alex)
        for i=1:4
            Visualizer(cur_best_solutions{i}, toy_nets{i}, threshold, transcendence)
        end
    end

    Toy_Solutions_Exp{sample} = best_solutions;

end

%     Alex: (i) run this N times (20ish) to get a distribution ofsolutions/time-to-best solution (num calls).
%     (ii) Compare this to a brute-force found optimal solution (~a few hrs
%     brute force time for EACH net of size 100 choose 3

% Toy_Mat is a Matrix wit four columns (fitness, # calls, trans value, network type
% 1 2 3 and 4 (lattice, star, chain and ER respectively))

%writematrix(ToyMat, 'toyNWdataFinal.csv')
toc
%csvwrite('toyNWdata1.csv',ToyMat)
%% (1c) Figure 1: Diagram of fitness calculation
%TODO
%probably outside of Matlab, Blake will do this


%NOTE: if we are to brute force search any of th

%% PART 2: ERDOS RENYI RANDOM GRAPH
%% (2a) create ER random graphs
%seed=27; %30, 32, 34, 50 for n=100 w/ 110,.025
%[ER_G, n, m]=ErdosRenyi(110, .02, seed);
%[ER_G, n, m]=ErdosRenyi(150,.012,randsample(1000,1));
[ER_G, n, m]=ErdosRenyi(53,.05,6);


%% (2b) run GA on the ER graph
rng('shuffle')
%initialize
P = 200; % GA population size
N = size(ER_G,1);
nGen=100;
clear V;
global V
V = 3; % # vaccines
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

% set options for ga toolbox 
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
    'PlotFcn',{@gaplotbestf},...
    'OutputFcn',{@outputFcn});

ER_Mat = zeros(N_reps, 2); % final ER matrix

for n = 1:N_reps

    % run GA
    [x,fval, exitFlag, Output] = ga(@(x) SpreadingFitnessFcnCompSize(x, ER_G, threshold, transcendence), V, vaccineOpts);
    record = outputFcn();
    fcncalls_to_best_solution = P * (record(end).LastImprovement + 1);
    clear outputFcn

    ER_Mat(n, 1) = fval;
    %ER_Mat(n, 2) = Output.funccount;
    ER_Mat(n, 2) = fcncalls_to_best_solution;
    
% Visualizer(x, ER_G, threshold, transcendence)

end

% TODO: run X times to get a distribution of fitnessess (& time to best)
% completed: ER_Mat is a matrix where the first column are the fitness
% values and second are number 0f function calls. Every row is a replicate
% run

%% (2b - part 2) Compare GA ER results to brute force search
brute_vaccs = nchoosek(1:N, V);
 
fval = zeros(size(brute_vaccs,1),1);
tic
for i = 1:size(brute_vaccs,1)
    fval(i) = SpreadingFitnessFcnCompSize(brute_vaccs(i, :), ER_G, ...
        threshold, transcendence);
end
toc

minimum = min(fval);
equal_best_solns = find(fval == minimum);

best_vacc = brute_vaccs(equal_best_solns, :);

% visually compare the GA and brute force true best solution
Visualizer(best_vacc(1, :), ER_G, threshold, transcendence);
saveas(gcf,'bruteForceSoln.fig');
Visualizer(x, ER_G, threshold, transcendence);
saveas(gcf,'gaSoln.fig');

%% create figure comparing true best solution and ga solution
% Load saved figures
b=hgload('bruteForceSoln.fig');
g=hgload('gaSoln.fig');

% Prepare subplots
figure
h(1) = subtightplot(1,2,1,[0.05,0.01], 0.075, 0.05);
set(gca,'XColor', 'none','YColor','none')
h(2) = subtightplot(1,2,2,[0.05,0.01], 0.075, 0.05);
set(gca,'XColor', 'none','YColor','none')

% Paste figures on the subplots
copyobj(allchild(get(b,'CurrentAxes')),h(1));
copyobj(allchild(get(g,'CurrentAxes')),h(2));

% Add legends
l(1)=title(h(1),'\fontsize{24}Global Optimum');
l(2)=title(h(2),'\fontsize{24}Best Solution Found');

saveas(gcf,'compareBruteGA.fig');

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

% completed by Thomas below

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
l(4)=title(h(4),'\fontsize{36}Erd?s?Rï¿½nyi');

%% PART 3: RUN ON FLU NETS OF VARIOUS SIZES FOR ~3 DIFFERENT TRANSCENDENCE VALUES
%% (3a) import flu networks
load('H3N2_flu_net_components.mat');
flu_nets=adjmats; % ordered from largest to smallest
clear adjmats

%% (3a pt2) brute force search flu networks
fvals={};
ct=1;
for net_idx=9:-1:5 %search the smallest 4 networks
    disp(net_idx)
    %get network
    cur_net=flu_nets{net_idx};
    
    N=size(cur_net,1);
    V=3; %3 vaccine case
    threshold=0.5;
    transcendence=1;
    
    %get all possible vaccination combinations
    brute_vaccs = nchoosek(1:N, V);
    %brute_vaccs = brute_vaccs(1:200,:);
    fval = zeros(size(brute_vaccs,1),1); %init fitness of each

    tic
    for i = 1:size(brute_vaccs,1)
        if mod(i,1000)==0
            disp(i)
        end
        fval(i) = SpreadingFitnessFcnCompSize(brute_vaccs(i, :), cur_net, ...
            threshold, transcendence);
    end
    toc
    fvals{ct}=fval;
    ct=ct+1;
end

brute_force_real_3vac=fvals;
save('bruteForceReal3vac.mat','brute_force_real_3vac')

 %% (3a pt3) visualize brute force searches
V=3;
threshold=.5;
transcendence=1;

load('bruteForceReal3vac.mat')
best_vaccs={};
for net_idx=1:length(brute_force_real_3vac)

    cur_fval=brute_force_real_3vac{net_idx};
    cur_minimum = min(cur_fval);
    equal_best_solns = find(cur_fval == cur_minimum);
    
    N=size(flu_nets{10-net_idx},1);
    brute_vaccs=nchoosek(1:N, V);
    
    cur_best_vacc = brute_vaccs(equal_best_solns, :);
    best_vaccs{net_idx}=cur_best_vacc;
    Visualizer(cur_best_vacc(1, :), flu_nets{10-net_idx}, threshold, transcendence);
    saveas(gcf,strcat('bruteForceSoln3vacRealN',int2str(N),'.fig'));

end

%% (3b) run GA on flu nets

num_GA_runs=20;
P = 300; % GA population size
nGen=50;
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
'PlotFcn',{@gaplotbestf},...
'OutputFcn',@outputFcn);


counter = 0;
transList = [.5 1 1.5];
RealMat = zeros((num_GA_runs * length(transList) * length(flu_nets)), 4);

for transcendence=1:length(transList) % since we don't know what it is, try a few
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
            
            counter = counter + 1;
            
            % run GA
            [x, fval, exitFlag, Output] = ga(@(y) SpreadingFitnessFcnCompSize(y, flu_net, threshold, transcendence), V, vaccineOpts);
            record = outputFcn();
            fcncalls_to_best_solution = P * (record(end).LastImprovement + 1);
            clear outputFcn
            %Visualizer(x, flu_net, threshold, transcendence)
            RealMat(counter,1) = fval;
            RealMat(counter,2) = fcncalls_to_best_solution;%Output.funccount;
            
            RealMat(counter,3) = transList(transcendence);
            RealMat(counter,4) = N;
        end

    end
end

%Done: Extract data for X number of runs of the above 
% RealMat matrix first col fitness the seoncond is number of func calls
% third is trans value and the fourth is the netwrok 1:9 small:big - Alex


%% (3c) Figure 3: Expected # strains in outbreak (fitness) vs network size (N), by transcendence.
% uses results of 3b

% writing out and and making figures and conducting stat anlaysis in R
% because Matlab is shit at figures and as statisitcal package - Alex
%csvwrite('realNWdata3.csv',RealMat)


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
l(1)=title(h(1),'\fontsize{24}N=20');
l(2)=title(h(2),'\fontsize{24}N=43');
l(3)=title(h(3),'\fontsize{24}N=52');
l(4)=title(h(4),'\fontsize{24}N=81');
l(5)=title(h(5),'\fontsize{24}N=105');
l(6)=title(h(6),'\fontsize{24}N=223');
l(7)=title(h(7),'\fontsize{24}N=400');
l(8)=title(h(8),'\fontsize{24}N=791');
l(9)=title(h(9),'\fontsize{24}N=1430');

% saved the example real nets into a 3 by 3 figure - Alex



%%  PART 4: RUN ON GROWING NETWORK FLU
%% (4a) import network

%import nodes
filename = 'TemporalNodes.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
temporalNodes = [dataArray{1:end-1}];
clearvars filename delimiter startRow formatSpec fileID dataArray ans;


%import edges
filename = 'TemporalEdges.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
temporalEdges = [dataArray{1:end-1}];
temporalEdges = sortrows(temporalEdges,3);
clearvars filename delimiter startRow formatSpec fileID dataArray ans;


%% (4b) Grow network and evaluate fitness at each step
nodeTimes=temporalNodes(:,2);
[minT maxT]=bounds(nodeTimes);

increment=50; %# days between evaluations.
ns=[];

days=sort([4772 minT:increment:maxT]);%from day of first observed strain to last
x=[]; %initialize solution
fvals_ga=[]; %initialize fitnesses
fvals_rand=[]; %initialize fitnesses
novacc=[];
%make full component
nodes=array2table((temporalNodes(:,1)));
G_full=graph();
G_full=addnode(G_full, nodes);
source=temporalEdges(:,1);
target=temporalEdges(:,2);
G_full=G_full.addedge(source,target);
G_full=full(adjacency(G_full));
ct=1;
population=[];
for day=days 
    if day>=4772 %if network is 1/2 full size by # nodes
        rng(4)
        %find edges in existence at this time
        idxs=temporalEdges(:,3)<=day;

        %make network
        source=temporalEdges(idxs,1);
        target=temporalEdges(idxs,2);
        G=graph(source, target);
        G=full(adjacency(G));
        [bins,binsizes] = conncomp(graph(G));

        %run GA on first half of net
        threshold = .5;
        transcendence=1;
        
        
        if day==4772
            ga_reps=10; %number of reps for the GA
            P = 200; % GA population size
            nGen=20;
            global V
            V = 4; % # vaccines
            mutProb = 1/V; % probabilty of mutation
            
            ga_solutions=zeros(ga_reps,V);
            
            for k = 1:ga_reps
            
                N = size(G,1);
                population=zeros(P,V);
                for j=1:P
                    population(j,:)=randsample(1:N,V);
                end
                
                
                vaccineOpts = gaoptimset(...
                    'PopulationType', 'doubleVector', ...
                    'PopulationSize',P,...
                    'InitialPopulation',population,...
                    'CrossoverFcn', @crossoversinglepoint, ...
                    'CrossoverFraction', 0.5, ...
                    'SelectionFcn',{@selectiontournament,2}, ... {@selectionstochunif}
                    'Vectorized','on',...
                    'Generations', nGen, ...
                    'FitnessLimit', 0, ...
                    'MutationFcn', {@randomResetMutation, mutProb},...
                    'PlotFcn',{@gaplotbestf},...
                    'OutputFcn',@outputFcn);
                [ga_solutions(k,:), fval(k), exitFlag, Output] = ga(@(y) SpreadingFitnessFcnCompSize(y, G, threshold, transcendence), V, vaccineOpts);
            
            end
            fvals_ga(1,:) = fval;
            
            % find the top 5% of random vaccination methods, testing as
            % many as population size * generations to best GA solution.
            record = outputFcn();
            test_number = P * (record(end).LastImprovement + 1);
            disp((record(end).LastImprovement + 1))
            clear outputFcn
            
            randvacc=zeros(test_number,V);
            for j=1:test_number
                randvacc(j,:)=randsample(1:N,V);
            end
            best_fitness=ones(1,round(test_number*.05));
            best_solutions=zeros(round(test_number*.05),V);
            for j=1:test_number
                curr_vacc=randvacc(j,:);
                curr_fit=SpreadingFitnessFcnCompSizeGrow(curr_vacc,G,G_full,threshold,transcendence);
                curr_worst=max(best_fitness);
                if curr_fit<curr_worst
                    idx=find(best_fitness==curr_worst,1);
                    best_fitness(idx)=curr_fit;
                    best_solutions(idx,:)=curr_vacc;
                end
            end
            fvals_rand=best_fitness;
        else
            %give it knowledge of the full network for fitness calcs,
            %negligible difference.
             fvals_ga_new=[];
            for j=1:size(ga_solutions,1)
                curr_vacc=ga_solutions(j,:);
                fvals_ga_new=[fvals_ga_new SpreadingFitnessFcnCompSizeGrow(curr_vacc,G,G_full,threshold,transcendence)];
            end
            fvals_ga=[fvals_ga;fvals_ga_new];           
          
            
            fvals_rand_new=[];
            for j=1:length(best_solutions)
                curr_vacc=best_solutions(j,:);
                fvals_rand_new=[fvals_rand_new SpreadingFitnessFcnCompSizeGrow(curr_vacc,G,G_full,threshold,transcendence)];
            end
            fvals_rand=[fvals_rand;fvals_rand_new];
        end
        ct=ct+1;
        
        [~, binsize]=conncomp(graph(G));
        disp(binsize)
        novacc=[novacc sum(binsize.^2)/(size(G,1)^2)];
        
    end
end

plot_days=days(days>=4772)-4772;
figure
hold on
plot(plot_days,fvals_rand,'-','linewidth',10,'color',[.5 .5 .5 .2])
plot(plot_days,fvals_ga,'-','linewidth',2)
plot(plot_days,novacc,':','linewidth',3)
hold off 




%% Test random vaccination strategies for toy NWs
randSize = 1000;
threshold = .5;
transcendenceList = [2 2.5 3];

%nw = toy_nets{1};

global V
V = 3; % # vaccines




randToyMat = zeros(length(toy_nets)*length(transcendenceList), 2+randSize);

counter=0; % intialize counter to index matrix


for transcend_idx=1:length(transcendenceList)
    transcendence=transcendenceList(transcend_idx);
    
        for i=1:length(toy_nets)

        %get network and N
        nw=toy_nets{i};
        N=size(nw,1);

        counter=counter+1; % counter to fill matrix
        
        % pick solutions 
        random_vaccs=zeros(randSize,V);
        for j=1:randSize
            random_vaccs(j,:)=randsample(1:N,V);
        end

        % evaluate solutions
        trials=SpreadingFitnessFcnCompSize(random_vaccs, nw, threshold, transcendence);

        randToyMat(counter,1) = i;
        randToyMat(counter,2) = transcendence;
        randToyMat(counter,3:3+randSize-1) = trials;
        
        end        

end

writematrix(randToyMat, 'randToyMat.csv')




%% Test random vaccination strategies for Real NWs (3 vac)
randSize = 1000;
threshold = .5;
transcendenceList = [1 2 3];


global V
V = 3; % # vaccines


randRealMat3 = zeros(length(flu_nets)*length(transcendenceList), 2+randSize);

counter=0; % intialize counter to index matrix


for transcend_idx=1:length(transcendenceList)
    transcendence=transcendenceList(transcend_idx);
    
        for i=1:length(flu_nets)

        %get network and N
        nw=flu_nets{i};
        N=size(nw,1);

        counter=counter+1; % counter to fill matrix
        
        % pick solutions 
        random_vaccs=zeros(randSize,V);
        for j=1:randSize
            random_vaccs(j,:)=randsample(1:N,V);
        end

        % evaluate solutions
        trials=SpreadingFitnessFcnCompSize(random_vaccs, nw, threshold, transcendence);

        randRealMat3(counter,1) = i;
        randRealMat3(counter,2) = transcendence;
        randRealMat3(counter,3:3+randSize-1) = trials;
        
        end        

end

writematrix(randRealMat3, 'randRealMat3.csv')


%% Test random vaccination strategies for Real NWs (4 vac)
randSize = 1000;
threshold = .5;
transcendenceList = [1 2 3];



global V
V = 4; % # vaccines




randRealMat4 = zeros(length(toy_nets)*length(transcendenceList), 2+randSize);

counter=0; % intialize counter to index matrix


for transcend_idx=1:length(transcendenceList)
    transcendence=transcendenceList(transcend_idx);
    
        for i=1:length(flu_nets)

        %get network and N
        nw=flu_nets{i};
        N=size(nw,1);

        counter=counter+1; % counter to fill matrix
        
        % pick solutions 
        random_vaccs=zeros(randSize,V);
        for j=1:randSize
            random_vaccs(j,:)=randsample(1:N,V);
        end

        % evaluate solutions
        trials=SpreadingFitnessFcnCompSize(random_vaccs, nw, threshold, transcendence);

        randRealMat4(counter,1) = i;
        randRealMat4(counter,2) = transcendence;
        randRealMat4(counter,3:3+randSize-1) = trials;
        
        end        

end

writematrix(randRealMat4, 'randRealMat4.csv')






