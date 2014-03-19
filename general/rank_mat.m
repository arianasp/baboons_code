function [ new_mat, ranks ] = rank_mat( mat, test_statistic, rank_direc, rank_dim)
%Puts the row and columns of an NxN matrix in rank order. Can do this 
%ranking by the values of the test statistic computed for each row 
%(default) or column. Can also do it in ascending order (default) or 
%descending order.
%INPUTS:
%   mat: [NxN matrix] of values to rank
%   test_statistic: [string] what test statistic to use as the value to 
%       rank by (can be 'mean','median','sum','min','max','std')
%   rank_direc: [string] specifying whether to rank in ascending order
%       or descending order ('ascend' or 'descend')
%   rank_dim: [string] specifying whether to rank by row or by column
%       ('row' or 'col'). If by row, the function treats each row as an
%       individual, takes the test statistic of each row, and does the
%       ranking based on this. Vice versa if by column
%OUTPUTS:
%   new_mat: reordered matrix
%   ranks: for each row (or column) in the original matrix, what is its
%       rank in the reordered matrix?
%DEFAULTS:
%   test_statistic = 'mean'
%   rank_direc = 'ascend'
%   rank_dim = 'row'
%MISSING DATA:
%   the function can handle rows / cols of NaNs - if using 'ascend', they 
%   are given the maximum rank (N). if using 'descend', they are given the 
%   minimum rank (1), consistent with what matlab's "sort" function 
%   returns. In other words, NaN is treated as the largest number in both 
%   cases. If there are multiple NaNs, their order is determined by their 
%   order in the original matrix (again consistent with "sort")

%matrix dimensions
N = size(mat,1);
M = size(mat,2);
if N~=M
    error('input matrix must be NxN')
end

%set defaults
if ~exist('test_statistic')
    test_statistic = 'mean';
end
if ~exist('rank_direc')
    rank_direc = 'ascend';
end
if ~exist('rank_dim')
    rank_dim = 'row';
end

%get the dimension to take the test statistic across
if strcmp(rank_dim,'row')
    dim = 2;
elseif strcmp(rank_dim,'col')
    dim = 1;
else
    error('must specify rank_dim as either "row" or "col"')
end

%get the values to rank on (based on the selected test statistic)
if strcmp(test_statistic,'mean')
    stats = nanmean(mat,dim);
elseif strcmp(test_statistic,'median')
    stats = nanmedian(mat,dim);
elseif strcmp(test_statistic,'sum')
    stats = nansum(mat,dim);
elseif strcmp(test_statistic,'min')
    stats = nanmin(mat,[],dim);
elseif strcmp(test_statistic,'max')
    stats = nanmax(mat,[],dim);
elseif strcmp(test_statistic,'std')
    stats = nanstd(mat,dim);
else
    error('unknown test statistic specified')
end

%sort ascending or descending to get ranks
if strcmp(rank_direc,'ascend') || strcmp(rank_direc,'descend')
    [vals order] = sort(stats,rank_direc);
else
    error('must specify rank_direc as either "ascend" or "descend"');
end

%put the matrix in the order specified by the ranks
new_mat = mat(order,order);

%make a vector 'ranks' which gives the rank of each row and col from the
%original matrix in the new matrix
ranks = order;
ranks(ranks) = 1:N;




end

