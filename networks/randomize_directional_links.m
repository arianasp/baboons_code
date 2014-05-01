function [ adj_mat ] = randomize_directional_links( adj_mat, missing_or_existing )
%For use with linearity_test
%Takes in an adjacency matrix with NaN's where missing links are. If
%missing_or_existing == 'missing', randomly puts in adj_mat(i,j) = 0 and 
%adj_mat(j,i) = 1 or vice versa in places where 
%adj_mat(i,j) = adj_mat(j,i) = NaN. If missing_or_existing == 'existing",
%does the same, but for non-NaN entries only.
%INPUTS:
%   adj_mat: [N x N matrix] adjacency matrix
%   missing_or_existing: [string] specifying whether to randomize
%       directional links between existing data (non-NaN entries) or
%       missing data (NaN entries)
%OUTPUTS;
%   adj_mat: [N x N matrix] randomized adjacency matrix

%find row and column indices of nans (or non-nans) in matrix, depending on
%whether you want to randomize existing links or misisng links
if strcmp(missing_or_existing,'missing')
    [row col] = find(isnan(adj_mat));
elseif strcmp(missing_or_existing,'existing')
    [row col] = find(not(isnan(adj_mat)));
else
    error('must specify missing or existing links')
end

%randomize missing links
for idx = 1:size(row,1)
    i = row(idx);
    j = col(idx);
    
    if i==j
        adj_mat(i,j)=0;
    else
        if rand < 0.5
            adj_mat(i,j) = 1;
            adj_mat(j,i) = 0;
        else
            adj_mat(i,j) = 0;
            adj_mat(j,i) = 1;
        end
    end
end
    


end

