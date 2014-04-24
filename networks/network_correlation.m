function [ R ] = network_correlation( A1, A2, corr_type, directed )
%Computes the correlation between two networks, described by adjacency
%matrix A1 and A2, using a correlation of type corr_type. If the networks
%are directed (directed == 1), the correlation is over all entries in the 
%matrices (excluding the diagonal). If they are not (directed == 0), then 
%only the upper triangle is used (to avoid double counting).
%For now, only the Pearson correlation is supported, but I may be adding
%more functionality later.
%INPUTS:
%   A1: [N x N matrix] of network connections 
%   A2: [N x N matrix] of network connections
%   corr_type: [string] type of correlation to use. Right now, supported
%   types are: 'Pearson'
%   directed: [bool] whether the network is directed (1) or undirected (0)
%OUTPUTS:
%   R: value of the correlation

N = size(A1,1);

if not(size(A1,2) == N && size(A2,1) == N && size(A2,2) == N)
    error('input matrices must be NxN and must be the same size as one another')
end

%get matrix elemenets in vector form
%if directed, use all terms except diagonal
%if undirected, use only upper triangular portion
uppertri = triu(ones(N),1);
lowertri = tril(ones(N),1);
if directed
    edges1 = [A1(find(uppertri)); A1(find(lowertri))];
    edges2 = [A2(find(uppertri)); A2(find(lowertri))];
else
    edges1 = A1(find(uppertri));
    edges2 = A2(find(uppertri));
end

if strcmp(corr_type,'Pearson') || strcmp(corr_type,'pearson')
    R = corrcoef(edges1,edges2,'rows','complete');
    R = R(1,2);
end
    
    


end

