%% create a test graph

N_prelim=40;
% set network density: 
p=.07;
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
immunity(n1:n1+2)=0; %partially vaccinate 3 nodes

adjmat=full(G.adjacency);
figure; plot(G)
%% create a test immunity

immunity = zeros(1,N);

%EDIT THE FOLLOWING LINE: set n1=1,2,3,4, up to N-2 etc
n1=1; %pick node indices n1, n1+1, and n1+2 to vaccinate 50% of the population for
immunity(n1:n1+2)=1; %partially vaccinate 3 nodes

disp(immunity)

%%
tic
[Ifinal convergence] = SIRSfitnessfcn(1000, 1/2,1/3.5, 1/1000,...
    2, 1/250, adjmat, immunity, 100000,10^-11);
toc
disp(['p(I)=' num2str(Ifinal,'%0.3g')])
if convergence==1
    disp('converged? YES')
else
    disp('converged? NO')
end