function [ data ] = mutual_information( X, Y, mins, maxes, nbins )
%Computes the mutual information between X and Y, which can
%be either Nx1 column vectors (i.e. one-dimensional data) or Nx2 matrices
%(i.e. two-dimensional data). Note that higher-dimensional data are not
%currently supported. This function takes in continuous data and first bins
%it in equally-sized bins according to nbins, which is either a single
%number (for one-dimensional data) or a 1x2 vector (for two-dimensional 
%data) which gives the number of bins to use for each dimension separately.
%INPUTS:
%   X: [Nx1 or Nx2 matrix] of 1-D or 2-D data
%   Y: [Nx1 or Nx2 matrix] of 1-D or 2-D data
%   mins: [number or 1x2 vector] giving the minimum value to include for
%       each dimension (matrix column)
%   maxes: [number or 1x2 vector] giving the maximum value to include for
%       each dimension (matrix column)
%   nbins: [number or 1x2 vector] giving the number of bins to use for each
%       dimension (matrix column)
%OUTPUTS:
%   data: [struct] containing the following fields
%       MI: [number] the mutual information between the two sets of data
%       H_X: entropy of the distribution determined by X
%       H_Y: entropy of the distribution determined by Y
%       n_times: number of non-NaN times that went into the computation

%Note: Deals with NaNs by removing rows of both X and Y where either X or Y
%contained a NaN in that row.

%Note: Removes data that is out of range (not within the range min to max
%in one or the other (or both) dimensions).

%make sure inputs are in the correct format
N = size(X,1);
D = size(X,2);

if size(Y,1) ~= N || size(Y,2) ~=D
    error('matrices X and Y must be the same size')
end

if D ~= 1 && D ~=2
    error('X and Y must be Nx1 or Nx2 matrices')
end

if length(nbins) ~= D || length(mins) ~= D || length(maxes) ~=D
    error('nbins, mins, and/or maxes vector does not correspond to the dimensionality of the X and Y matrices')
end

if N <= 2
    warning('2 or fewer data points in X and Y matrices')
end

%remove NaNs
nansx1 = find(isnan(X(:,1)));
nansy1 = find(isnan(Y(:,1)));
nans = union(nansx1,nansy1);
if D == 2
    nansx2 = find(isnan(X(:,2)));
    nansy2 = find(isnan(Y(:,2)));
    nans = union(nans,nansx2);
    nans = union(nans,nansy2);
end
X(nans,:) = [];
Y(nans,:) = [];

%bins for dimension 1
x1 = transpose(linspace(mins(1), maxes(1),nbins(1)));

%interpolaed values for x and y in dimension 1
xr1 = interp1(x1,1:numel(x1),X(:,1),'nearest');
yr1 = interp1(x1,1:numel(x1),Y(:,1),'nearest');

if D == 2
    %bins for dimension 2
    x2 = transpose(linspace(mins(2), maxes(2),nbins(2)));
    %interpolated values for x and y in dimension 2
    xr2 = interp1(x2,1:numel(x2),X(:,2),'nearest');
    yr2 = interp1(x2,1:numel(x2),Y(:,2),'nearest');
end

%remove out of range data
nans = find(isnan(xr1));
nans = union(nans,find(isnan(xr2)));
nans = union(nans,find(isnan(yr1)));
nans = union(nans,find(isnan(yr2)));
xr1(nans) = [];
xr2(nans) = [];
yr1(nans) = [];
yr2(nans) = [];

n_times = length(xr1);

%construct distributions
if D == 2
    joint = [xr1, yr1, xr2, yr2];
    sing_x = [xr1, xr2];
    sing_y = [yr1, yr2];
    joint = accumarray(joint,1,[nbins(1) nbins(1) nbins(2) nbins(2)]);
    sing_x = accumarray(sing_x,1,[nbins(1) nbins(2)]);
    sing_y = accumarray(sing_y,1,[nbins(1) nbins(2)]);
    joint = joint / sum(sum(sum(sum(joint))));
    sing_x = sing_x / sum(sum(sing_x));
    sing_y = sing_y / sum(sum(sing_y));
elseif D == 1
    joint = [xr1,yr1];
    sing_x = xr1;
    sing_y = yr1;
    joint = accumarray(joint,1,[nbins nbins]);
    sing_x = histc(sing_x,x1);
    sing_y = histc(sing_y,x1);
    joint = joint / sum(sum(joint));
    sing_x = sing_x / sum(sing_x);
    sing_y = sing_y / sum(sing_y);
end

%get log values and substitute 0 for log(0) so that 0*log(0) = 0
log_joint = log(joint);
log_sing_x = log(sing_x);
log_sing_y = log(sing_y);
log_joint(find(joint == 0)) = 0;
log_sing_x(find(sing_x == 0)) = 0;
log_sing_y(find(sing_y == 0)) = 0;

%compute entropy
if D == 1
    H_joint = -sum(sum(joint.*log_joint));
    H_x = -sum(sing_x.*log_sing_x);
    H_y = -sum(sing_y.*log_sing_y);
elseif D == 2
    H_joint = -sum(sum(sum(sum(joint.*log_joint))));
    H_x = -sum(sum(sing_x.*log_sing_x));
    H_y = -sum(sum(sing_y.*log_sing_y));
end

%compute mutual information
MI = H_x + H_y - H_joint;

data.MI = MI;
data.H_joint = H_joint;
data.H_X = H_x;
data.H_Y = H_y;
data.n_times = n_times;






    




end

