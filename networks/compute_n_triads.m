function [ d ] = compute_n_triads( adj_mat )
%Computes the number of cyclical triads in an adjacency matrix adj_mat.
%linearity_test calls this function.
%this function is based on de vries 1995
%adj_mat must be of the form adj_mat(i,j) = 1 and adj_mat(j,i) = 0 if i
%dominates j (in other words, if adj_mat(i,j) = 1 then adj_mat(j,i) must be
%0, except on the diagonals which are 0 (adj_mat(i,i) = 0)
%INPUTS:
%   adj_mat: [NxN matrix] adjacency matrix
%OUTPUTS;
%   d: number of cyclical triads in the network specified by adj_mat


if not(size(adj_mat,1) == size(adj_mat,2))
    error('requires a square matrix')
end

%number of nodes
N = size(adj_mat,1);

%sums along the rows
Si = sum(adj_mat,2);

%compute the number of triads
d = N*(N-1)*(2*N-1)/12 - 1/2*sum(Si.^2);


end

