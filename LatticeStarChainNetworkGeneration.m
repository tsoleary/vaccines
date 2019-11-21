% GENERATE LATTICE, STAR, AND CHAIN NETWORKS 
% Nx: lattice x dimension 
% Ny: lattice y dimension
% N=Nx*Ny, so N is constant among all graph types
Nx=8;
Ny=8;
N=Nx*Ny;

% make lattice
lattice_adjmat=zeros(Nx,Ny);
for i=1:numel(lattice_adjmat)
    if mod((i-1),Nx)~=0
        lattice_adjmat(i,i-1)=1;
        lattice_adjmat(i-1,i)=1;
    end
    if mod(i,Nx)~=0
        lattice_adjmat(i,i+1)=1;
        lattice_adjmat(i+1,i)=1;
    end
    if i>Nx
        lattice_adjmat(i,i-Nx)=1;
        lattice_adjmat(i-Nx,i)=1;
    end
    if i<(Nx*Ny-Nx)
        lattice_adjmat(i,i+Nx)=1;
        lattice_adjmat(i,i+Nx)=1;
    end
end

% make star
star_adjmat=zeros(N);
star_adjmat(2:end,1)=ones(N-1,1);
star_adjmat(1,2:end)=ones(N-1,1);

% make chain
chain_adjmat=zeros(N);
for i=2:N-1
    chain_adjmat(i,i-1)=1;
    chain_adjmat(i,i+1)=1;
    chain_adjmat(i+1,i)=1;
    chain_adjmat(i-1,i)=1;
end

lattice_adjmat;
star_adjmat;
chain_adjmat;

%figure; plot(graph(lattice_adjmat),'layout','force');
%figure; plot(graph(star_adjmat),'layout','force');
%figure; plot(graph(chain_adjmat),'layout','force');
