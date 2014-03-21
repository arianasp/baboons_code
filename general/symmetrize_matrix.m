function [ symm_mat ] = symmetrize_matrix( mat )
%Takes in a matrix of counts of directed interactions, and returns a matrix
%of counts of all interactions. 
%So symm_mat(i,j) = symm_mat(j,i) = mat(i,j) + mat(j,i)
%INPUTS:
%   mat: [NxN matrix] 
%OUTPUTS:
%   symm_mat: [NxN matrix] defined as above

N = size(mat,1);
if not(size(mat,2)==N)
    error('must input a symmetric matrix')
end

symm_mat = nan(N,N);

for i = 1:N
    for j = 1:N
        symm_mat(i,j) = mat(i,j) + mat(j,i);
        symm_mat(j,i) = mat(i,j) + mat(j,i);
    end
end


end

