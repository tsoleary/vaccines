%%% README

% This file:
% (1) generates an Erdos-Renyi random graph of specified size and density
% (2) included immunization by reducing transmission rates at certain nodes
% (3) runs the graph through an SIRS model 




%%% SET PARAMETERS

rng(3)

% set network size:
N_prelim=20;

% set network density: 
p=.2;
%k_avg=2;
%p=k_avg * 2/(N*(N-1)); 

% immunity set later


%%% create ER random graph

% network density


%create fully connected edgelist
full_edgelist = [];
for i=1:N_prelim
    for j=i+1:N_prelim
        full_edgelist=[full_edgelist;i j];
    end
end

% find # edges to remove
edgelist=full_edgelist;
n_end=round(p * N_prelim * (N_prelim-1)/2);
n_remove=size(full_edgelist,1)-n_end;

%random shuffle edges
edgelist_shuffle = edgelist(randperm(size(edgelist, 1)), :);

%remove n edges
size(edgelist)
edgelist=edgelist_shuffle(1:end-n_remove,:);
size(edgelist)

%insert edges into adjacency matrix
adjmat=zeros(N_prelim);
for k=1:size(edgelist,1)
    x=edgelist(k,1);
    y=edgelist(k,2);
    adjmat(x,y)=1;
    adjmat(y,x)=1;
end

% create graph
G_prelim=graph(adjmat);

% restrict to giant component, to ensure connectedness
[bin,binsize] = conncomp(G_prelim);
G = subgraph(G_prelim, binsize(bin) == max(binsize));

% update network size
N=G.numnodes;

% set immunity: 
immunity = zeros(1,N);

%EDIT THE FOLLOWING LINE: set n1=1,2,3,4, up to N-2 etc
n1=1; %pick node indices n1, n1+1, and n1+2 to vaccinate 50% of the population for
immunity(n1:n1+2)=.5; %partially vaccinate 3 nodes

figure
subplot(1,2,1)
h=plot(G,'MarkerSize',15,'NodeColor','k','LineWidth',3,'Layout','force','EdgeColor','k');
highlight(h,n1:n1+2,'NodeColor',[.5 0 0])
for i=n1:n1+2
    nbors=neighbors(G,i);
    highlight(h,i,nbors,'EdgeColor',[.5 0 0],'LineWidth',2)

end
title('Strain network')



%%% set SIRS parameters

delta=2;% immunity transcendence: 0 means none, recommended values: {1,2, 0 up to 10}
tmax=400; % duration of simulation: set such that infection curves settle
alpha=1/2; %transmission rate (should be greater than recovery rate)
beta=1/3.5; %recovery rate
gamma=1/1000; %immunity loss rate 
epsilon=1/250; %mutation (influences spread throughout network. If infections tend to ->0 without immunization, increase this or play with network density.)
tspan=1:tmax;

S0=100000; % human population size
R0=zeros(1,N); % nobody starts with specific immunity (immunization implemented w/infection rate alteration)




%%% create outbreak at each node
SIRS_results={}; % store infection counts for each outbreak
totals=zeros(tmax,N);
steady_state_means=zeros(1,N); % store FINAL infection counts

tic
for j=1:N % for each strain
    disp(['node ' num2str(j)]) 
    distance_mat=distances(G); % shortest paths between nodes
    
    I0=zeros(1,N); % Nobody's infected with any strain, 
    I0(j)=S0/1000; % ... except for the outbreak strain.

    %actually runs it 
    [t,y]=ode113(@(t,y) SIRS(t,y, alpha, beta, gamma, delta, epsilon, distance_mat, immunity), tspan, [S0 I0 R0]);

    % retrieve infected counts (columns of y are [S I1 I2 ... In R1 R2 ... Rn], rows are time steps, values = # of population in that state at that time)
    n=floor((size(y,2)-1)/2);
    I_counts=y(:,2:n+1);
    SIRS_results{j}=I_counts;
    totals(:,j)=sum(I_counts,2);
    steady_state_means(j)=sum(I_counts(end,:));
    disp(sum(y(1,:)))
    disp(sum(y(round(tmax/2),:)))
    disp(sum(y(end,:)))
    
end
toc
disp(['N=' num2str(G.numnodes) ', M=' num2str(G.numedges)]);

clear h

%plot infections by node
subplot(1,2,2)
hold on
for i=1:N
    h(i)=plot(tspan,sum(SIRS_results{i},2)/S0,'color',[1 0 0 .5],'DisplayName','Outbreak by node');
end
h(N+1)=plot(tspan,mean(totals,2)/S0,'color','k','linewidth',3,'DisplayName','Mean');
hold off
legend(h(N:N+1));
xlabel('t')
ylabel('\% population infected ($$\sum_i I_i/N$$)','Interpreter', 'Latex')
% look at strain-specific infection counts for a handful of outbreaks 

%pick random strain outbreaks to look at
strain_idxs=datasample(1:N,5);
%
figure
clear h
for i=1:length(strain_idxs)
    strain_idx=strain_idxs(i);
    subplot(1,5,i)
    hold on
    h(1:N)=plot(tspan,SIRS_results{strain_idx},'color',[.5 .5 .5],'DisplayName','Other strains');
    h(N+1)=plot(tspan,SIRS_results{strain_idx}(:,strain_idx),'color',[.6 0 0],'linewidth',2,'DisplayName','Outbreak strain');
    legend(h(N:N+1))
    title({['Outbreak at strain ' num2str(strain_idx)] , 'Strain-specific infections'})
    if strain_idx>=n1 && strain_idx<=(n1+2)
        title({['Outbreak at strain ' num2str(strain_idx)] , 'Strain-specific infections','NOTE: VACCINATED OUTBREAK STRAIN'})
    end
    hold off

end


