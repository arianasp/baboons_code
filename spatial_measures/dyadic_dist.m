function [ dyad_dist ] = dyadic_dist( xs, ys, test_statistic )
%Computes the distance between every pair of individuals. Then returns data
%depending on the choice of test statistic.
%INPUTS:
%   xs: N x T matrix giving x positions of all individuals over time
%   ys: N x T matrix giving y positions of all individuals over time
%   test_statistic: indicates which test statistic of the resulting
%       distribution of dyadic distances to return. This can be 'mean',
%       'median', 'mean_per_ind' (to return the mean distance to all others
%       from each individual) or 'full' (to return all dyadic distances)
%       NOTE: 'mean' or 'median' will return the mean (or median) of the 
%       individual means (or medians), not the mean (or median) of all
%       values. This means that each individual is weighted equally.
%OUTPUTS:
%   dyad_dist: the specified statistic of the distribution. If
%       test_statistic = 'mean' or 'median', this will be a 1 X T vector,
%       if test_statistic = 'mean_per_ind', this will be a N x T array
%       if test_statistic = 'full', this will be a N x N x T array

n_inds = size(xs,1);
n_times = size(xs,2);


%defaults to test_statistic = 'mean'
if nargin == 2
    test_statistic = 'mean';
end

if strcmp(test_statistic,'mean_per_ind')
    dyad_dist = nan(n_inds,n_times);
elseif strcmp(test_statistic,'full')
    dyad_dist = nan(n_inds,n_inds,n_times);
else
    dyad_dist = nan(1,n_times);
end

for t = 1:n_times
    t
    currdat = [xs(:,t) ys(:,t)];
    dists = squareform(pdist(currdat,'euclidean'));
    dists(logical(eye(size(dists)))) = NaN; %NaNs on the diagonal
    if strcmp(test_statistic,'mean_per_ind')
        dists = nanmean(dists,2);
        dyad_dist(:,t) = dists;
    elseif strcmp(test_statistic,'full')
        dyad_dist(:,:,t) = dists;
    elseif strcmp(test_statistic,'mean')
        dyad_dist(t) = nanmean(nanmean(dists));
    elseif strcmp(test_statistic,'median')
        dyad_dist(t) = nanmedian(nanmedian(dists));
    end
end
    


end

