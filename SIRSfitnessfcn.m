function [Ifinal convergence] = SIRSfitnessfcn(tmax, alpha, beta, gamma,...
    delta, epsilon, contact, immunity, popsize, threshold)

tspan=1:tmax;
N=size(contact,1);

S0=popsize*.999; % human population size
R0=zeros(1,N); % nobody starts with specific immunity (immunization implemented w/infection rate alteration)

%0 means run longer, 1 means good
convergence=0;

%%% create outbreak at each node
SIRS_results={}; % store infection counts for each outbreak
totals=zeros(tmax,N);
steady_state_means=zeros(1,N); % store FINAL infection counts
near_steady_state_means=zeros(1,N); % store near-FINAL infection counts

tic
for j=1:1 % for each strain
    disp(['node ' num2str(j)]) 
    distance_mat=distances(graph(contact)); % shortest paths between nodes
    
    I0=ones(1,N)/N;% uniform outbreak. All strains become prevalent at 
                   % once. It appears that the steady state is equivalent 
                   % as single strain outbreaks.

    %actually runs it 
    [t,y]=ode113(@(t,y) SIRS(t,y, alpha, beta, gamma, delta, epsilon, distance_mat, immunity), tspan, [S0 I0 R0]);

    % retrieve infected counts (columns of y are [S I1 I2 ... In R1 R2 ... Rn], rows are time steps, values = # of population in that state at that time)
    n=floor((size(y,2)-1)/2);
    I_counts=y(:,2:n+1);
    SIRS_results{j}=I_counts;
    totals(:,j)=sum(I_counts,2);
    steady_state_means(j)=sum(I_counts(end,:));
    near_steady_state_means(j)=sum(I_counts(end-10,:));
    disp(sum(y(1,:)))
    disp(sum(y(round(tmax/2),:)))
    disp(sum(y(end,:))) 
    
    
end

if abs( mean(near_steady_state_means) - mean(steady_state_means))/ steady_state_means(j) < popsize*threshold
    convergence=1;
end

Ifinal=mean(steady_state_means)/popsize;
