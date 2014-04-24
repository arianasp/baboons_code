function [ A_perm ] = permute_edges( A, net_type )
%Permutes the edges in the network specified by the adjacency matrix A. The
%matrix can be directed, in which case all of the entries are shuffled,
%undirected in which the edges are shuffled by symmetry of the adjacency
%matrix is preserved, or antisymmetric (i.e. A(i,j) = -A(j,i)) in which 
%case the anti-symmetry is preserved. Any NaNs in the matrix will retain
%their original spots.
%INPUTS:
%   A: [N x N matrix] adjacency matrix, can be weighted or not, directed or
%       not
%   net_type: [string] specifying what type of network A is
%       'directed': directed network - all edges will be shuffled
%       'undirected' / 'symmetric' : undirected network - shuffling will
%           maintain undirected nature of edges (symmetric adjacency
%           matrix)
%       'anti-symmetric': anti-symmetric network where A(i,j) = -A(j,i).
%           This structure will be maintained in the permutation
%OUTPUTS:
%   A_perm: [N x N matrix] Permuted version of A

N = size(A,1);
if not(size(A,2)==N)
    error('A must be an NxN matrix')
end

switch net_type
    case 'directed'
        non_nans = find(not(isnan(A)));
        edges = A(non_nans);
        edges_shuff = edges(randperm(length(edges)));
        A(non_nans) = edges_shuff;
    case {'undirected','symmetric'}
        uppertri = triu(ones(N),1); 
        all_edges = A(find(uppertri));
        non_nans = find(not(isnan(all_edges)));
        edges = all_edges(non_nans);
        edges_shuff = edges(randperm(length(edges)));
        all_edges(non_nans) = edges_shuff;
        A = squareform(all_edges,'tomatrix');
        A(find(logical(eye(N)))) = NaN; 
    case {'anti-symmetric','antisymmetric'}
        uppertri = triu(ones(N),1); 
        lowertri = tril(ones(N),-1);
        all_edges = A(find(uppertri));
        non_nans = find(not(isnan(all_edges)));
        edges = all_edges(non_nans);
        edges_shuff = edges(randperm(length(edges)));
        all_edges(non_nans) = edges_shuff;
        A = squareform(all_edges,'tomatrix');
        A(find(logical(eye(N)))) = NaN;
        A(find(lowertri)) = - A(find(lowertri));
        
    otherwise
        error('net_type not recognized - possible net types are directed, undirected, and anti-symmetric')
end

A_perm = A;

end

