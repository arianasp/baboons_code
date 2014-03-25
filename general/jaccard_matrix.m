function [ J ] = jaccard_matrix( M )
%Compute the "jaccard matrix" from a given matrix, M.
%The Jaccard matrix (which is something I'm making up here, though I'm sure
%it exists in the literature already somewhere), is defined as
%J(i,j) = (M(i,j) + M(j,i)) / (sum(M(i,:)) + sum(M(:,j)))
%This is a way of normalizing to what fraction of the interactions that i
%and j have with anyone are interactions with each other.
%INPUTS:
%   M: [NxN matrix] probably representing the number of interactions 
%       between N individuals
%OUTPUTS:
%   J: [NxN matrix] the "Jaccard matrix" as defined above
%MISSING DATA:
%   ignores NaNs

N = size(M,1);

if size(M,2) ~= N
    error('the entered matrix must be NxN')
end


J = nan(N,N);

for i = 1:N
    for j = 1:N
        J(i,j) = M(i,j) + M(j,i) / (nansum(M(i,:)) + nansum(M(:,j)));
    end
end


end

