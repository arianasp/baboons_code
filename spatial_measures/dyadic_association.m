function [ data ] = dyadic_association( X, Y, mins, maxes, nbins )
%Computes the dyadic association between X and Y, which can
%be either Nx1 column vectors (i.e. one-dimensional data) or Nx2 matrices
%(i.e. two-dimensional data). 
%The dyadic association in a given bin (i,j) is given by:
%	A(i,j) = log(P_AB(i,j) / (P_A(i,j)*P_B(i,j))), 
%		where P_AB(i,j) is the probability that both A and B are found in
%		bin (i,j), P_A(i,j) is the probability that individual A is found in
%		bin (i,j), and P_B(i,j) is the probability that individual B is found
%		in bin (i,j).
%The total dyadic association is given by the weighted sum:
%	A_tot = sum(A(i,j)*P_A(i,j)*P_B(i,j))
%This function takes in continuous data and first bins
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
%       A_ij: [vector length n_bins or n_bins(1)xn_bins(2) matrix] with
%			A(i,j) values for each bin (i,j), as defined above
%		A_tot: [number] A_tot as defined above
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
x1 = transpose(linspace(mins(1)-(maxes(1)-mins(1))/(nbins(1)-2)/2, maxes(1)+(maxes(1)-mins(1))/(nbins(1)-2)/2,nbins(1)));

%interpolaed values for x and y in dimension 1
xr1 = interp1(x1,0:(numel(x1)-1),X(:,1),'nearest');
yr1 = interp1(x1,0:(numel(x1)-1),Y(:,1),'nearest');

if D == 2
    %bins for dimension 2
    x2 = transpose(linspace(mins(2)-(maxes(2)-mins(2))/(nbins(2)-2)/2, maxes(2)+(maxes(2)-mins(2))/(nbins(2)-2)/2,nbins(2)));
    %interpolated values for x and y in dimension 2
    xr2 = interp1(x2,0:(numel(x2)-1),X(:,2),'nearest');
    yr2 = interp1(x2,0:(numel(x2)-1),Y(:,2),'nearest');
end

%remove out of range data
oor_x1 = union(find(X(:,1)<mins(1)),find(X(:,1)>maxes(1)));
oor_y1 = union(find(Y(:,1)<mins(1)),find(Y(:,1)>maxes(1)));
oor = union(oor_x1,oor_y1);

if D == 2
    oor_x2 = union(find(X(:,2)<=mins(2)),find(X(:,2)>=maxes(2)));
    oor_y2 = union(find(Y(:,2)<=mins(2)),find(Y(:,2)>=maxes(2)));
    oor2 = union(oor_x2,oor_y2);
    oor = union(oor,oor2);
end

xr1(oor) = [];
yr1(oor) = [];

if D == 2
    xr2(oor) = [];
    yr2(oor) = [];
end

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

%get probability distributions for computing association
P_X = sing_x;
P_Y = sing_y;

P_XY = nan(nbins(1),nbins(2));
for i = 1:nbins(1)
	for j = 1:nbins(2)
		P_XY(i,j) = joint(i,j,i,j);
	end
end

A_ij = log(P_XY) - log(P_X) - log(P_Y);
A_ij(find(P_XY == 0 .* P_X == 0 .* P_Y == 0)) = NaN;
A_tot = nansum(nansum(A_ij.*P_X.*P_Y));

data.A_ij = A_ij;
data.A_tot = A_tot;
data.P_XY = P_XY;
data.P_X = P_X;
data.P_Y = P_Y;
data.n_times = n_times;






    




end

