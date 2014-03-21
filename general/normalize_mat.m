function [ norm_mat ] = normalize_mat( mat, norm_dim, norm_stat )
%Normalize a matrix by dividing each entry in a given row by the sum of all
%of the entries in that row (or similarly for columns).
%INPUTS:
%   mat: [NxM matrix]
%   norm_dim: [string] dimension to normalize on 
%       ('rows'/'row' or 'cols'/'col')
%   norm_stat: what to normalize based on (what to divide each row/column
%       by - defaults to 'sum' but can also be 'mean', 'median', 
%       'max' or 'min')
%OUTPUTS:
%   norm_mat: normalize matrix,  as described above
%Note about NaNs: ignores NaNs in creating sums, but preserves existing NaN
%entries in the matrix

N = size(mat,1);
M = size(mat,2);

norm_mat = nan(N,M);

if strcmp(norm_dim,'row') || strcmp(norm_dim,'rows')
    dim = 2;
elseif strcmp(norm_dim,'col') || strcmp(norm_dim,'cols')
    dim = 1;
else
    error('unknown normalization dimension')
end

if strcmp(norm_stat,'sum')
    denoms = nansum(mat,dim);
elseif strcmp(norm_stat,'mean')
    denoms = nanmean(mat,dim);
elseif strcmp(norm_stat,'median')
    denoms = nanmedian(mat,dim);
elseif strcmp(norm_stat,'max')
    denoms = nanmax(mat,[],dim);
elseif strcmp(norm_stat,'min')
    denoms = nanmin(mat,[],dim);
else
    error('unknown normalization statistic')
end

if strcmp(norm_dim,'row') || strcmp(norm_dim,'rows')
    for i = 1:N
        norm_mat(i,:) = mat(i,:) / denoms(i);
    end
elseif strcmp(norm_dim,'col') || strcmp(norm_dim,'cols')
    for i = 1:N
        norm_mat(:,i) = mat(:,i) / denoms(i);
    end
end


end

