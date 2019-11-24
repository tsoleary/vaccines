function [G,n,m] = ErdosRenyi(n, p, random_seed)

rng(random_seed)

% create random adjacency matrix
ER_prelim=triu(sprand(n,n,p),1);
ER_prelim=spones(ER_prelim);
ER_prelim=ER_prelim+ER_prelim';
G=graph(ER_prelim);
% extract giant component
[bin,binsize]=conncomp(G);
idx=binsize(bin)==max(binsize); 
GC=subgraph(G, idx);

G= adjacency(GC);
n=GC.numnodes;
m=GC.numedges;

figure; plot(graph(G),'layout','force')