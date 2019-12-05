
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




%% (3b) run GA on flu nets to compare to brute force

num_GA_runs=20;
P = 300; % GA population size
nGen=50;
global V
V = 3; % # vaccines
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
transList = [1];
RealMat = zeros((num_GA_runs * length(transList) * length(flu_nets)), 4);

for transcendence=1:length(transList) % since we don't know what it is, try a few
    for i=6:9 %for flu networks of size ~20 up to ~1000

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
% third is trans val
writematrix(RealMat, 'brutForceGAfig.csv')




